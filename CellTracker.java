import java.util.ArrayList;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.Overlay;
import ij.gui.PolygonRoi;
import ij.gui.Roi;
import ij.plugin.ImageCalculator;
import ij.plugin.PlugIn;
import ij.plugin.RoiEnlarger;
import ij.plugin.frame.RoiManager;
import ij.process.ByteProcessor;
import ij.process.ImageProcessor;
import ij.Prefs; 
import ij.WindowManager; 
import ij.gui.Overlay; 
import ij.measure.ResultsTable;
// No need to import Spot if it is in the same folder

public class CellTracker implements PlugIn {
	@Override
	public void run(String arg) {
		
		ImagePlus original = IJ.getImage();
		ImagePlus aligned_stack_cells = null; 
		ImagePlus aligned_stack_wells = null; 
		
		// During development we directly work on the aligned images
		boolean align = false; 
		if (align) {
			// Meaning of the numbers in Duplicator:
			// firstC, lastC, firstZ, lastZ, firstT, lastT
			// For experiment 1, the slice 2 is better for the channel detection
			int bestZ_wells = 2;
			int bestZ_cells = 1; 
			int lastT = original.getNFrames();
			// Extract only the second slice in the image
			aligned_stack_wells = new ij.plugin.Duplicator().run(original, 1, 1, bestZ_wells, bestZ_wells, 1, lastT); 
			aligned_stack_cells = new ij.plugin.Duplicator().run(original, 1, 1, bestZ_cells, bestZ_cells, 1, lastT); 
			// Align the images across time
			aligned_stack_wells.show(); 
			IJ.run(aligned_stack_wells, "StackReg ", "transformation=Affine");
			IJ.run(aligned_stack_cells, "StackReg ", "transformation=Affine");
		} else {
			aligned_stack_wells = IJ.openImage("/Users/quentindevaud/Desktop/EPFL/Master/MA2/Bioimage informatics/Mini_project/images/p1_stabilised_cells.tif");;
			aligned_stack_cells = IJ.openImage("/Users/quentindevaud/Desktop/EPFL/Master/MA2/Bioimage informatics/Mini_project/images/p1_stabilised_cells.tif");;
			aligned_stack_wells = original;
			aligned_stack_cells = original;
			
		}
		//aligned_stack_cells.show();
		
		RoiManager rm = track_wells(aligned_stack_wells); 
		ImagePlus cells_only = remove_background(rm, aligned_stack_wells, aligned_stack_cells); 
		// TODO improve the tracking of the cells
		ArrayList<ArrayList<ArrayList<Double>>> cells_position = track_cells(cells_only); // Position of the cells in terms of numbers of frames, X and Y, number of cells per frame
		draw_wells(rm, cells_position); // Draw color if cell is inside well
		

		//ImagePlus imp_ch = imp.crop("slice=2, frames=1-270");
		//imp_ch.setTitle("test"); 
		//imp_ch.show(); 
		//IJ.selectWindow("test");
		//IJ.run("Duplicate...", "duplicate slices=2");
		//IJ.run(imp_ch, "Find Edges", "stack");
		//IJ.run("Threshold...");
		//IJ.setAutoThreshold(imp_ch, "Minimum dark");
		//IJ.run("Threshold...");
		//IJ.run(imp_ch, "Analyze Particles...", "size=150-Infinity circularity=0.00-0.13 show=[Overlay Masks] clear add stack");
		//imp_ch.show();
		
		
		/*
		 * double sigma = 5; double threshold = 10; ImagePlus imp = IJ.getImage(); int
		 * nt = imp.getNFrames(); int nx = imp.getWidth(); int ny = imp.getHeight();
		 * 
		 * ImagePlus dog = dog(imp, sigma); ArrayList<Spot> localmax[] = localMax(dog);
		 * ArrayList<Spot> spots[] = filter(dog, localmax, threshold);
		 * 
		 * // New parameters double d_max = Math.sqrt(ny*ny+nx*nx); double s_max =
		 * Math.pow(2, imp.getBitDepth()); double lambda = 0.5; double thres = 1e-2;
		 * boolean nearest = true;
		 * 
		 * 
		 * // Loop on the frames for (int t = 0; t < nt - 1; t++) { // Loop on the spots
		 * of the frame t for (Spot current : spots[t]) { // Loop on the spots of the
		 * next frame for (Spot next : spots[t+1]) { // Nearest neighbour linking: first
		 * store all the potential neighbours if(nearest) current.add_neighbour(next);
		 * // Simple linking: link only if the cost is below a threshold else if
		 * (current.distance(next, lambda, d_max, s_max) < thres) current.link(next); //
		 * Then get the nearest neighbour among all the potential neighbours if(nearest)
		 * current.get_nearest_neighbour(lambda, d_max, s_max, thres); } } } Overlay
		 * overlay = new Overlay(); draw(overlay, spots); imp.setOverlay(overlay);
		 */
		
		
	}
	
	
	public ArrayList<ArrayList<ArrayList<Double>>> track_cells(ImagePlus imp) {
			
			// Get a mask with only the cells
			imp.setDisplayRange(2502, 2876); // TODO: Need automatic computation of these values
			IJ.run(imp, "Apply LUT", "stack");
			IJ.setAutoThreshold(imp, "Default no-reset");
			Prefs.blackBackground = false;
			IJ.run(imp, "Convert to Mask", "method=Default background=Light calculate");
			ImagePlus imp2 = imp.duplicate();
			// Get the DOG 
			IJ.run(imp, "Gaussian Blur...", "sigma=7.05 stack"); // Not used
			IJ.run(imp2, "Gaussian Blur...", "sigma=5 stack");
			//ImagePlus imp1 = WindowManager.getImage("p1_stabilised_cells.tif");
			//ImagePlus imp2 = WindowManager.getImage("p1_stabilised_cells-1.tif");
			//ImagePlus imp3 = ImageCalculator.run(imp2, imp, "Subtract create stack");
			// Make binary and get only the edges of the cells
			//IJ.setAutoThreshold(imp3, "Default no-reset");
			//IJ.run(imp3, "Convert to Mask", "method=Default background=Light calculate");
			//IJ.run(imp3, "Analyze Particles...", "size=200-Infinity circularity=0.35-1.00 show=Overlay clear add stack");
			//imp3.show();
			int nt = imp2.getNFrames();
			ArrayList<ArrayList<ArrayList<Double>>> Values_RT = new ArrayList<ArrayList<ArrayList<Double>>>();
			for (int t = 0; t < nt-1; t++) { // loop over frames
				ArrayList<ArrayList<Double>> matrix = new ArrayList<ArrayList<Double>>();
				imp2.setSlice(t); // select slice
				IJ.run(imp2, "Find Maxima...", "prominence=10 output=List");
				ResultsTable RT1 = ResultsTable.getResultsTable();
				ArrayList<Double> Xs = new ArrayList<Double>();
				ArrayList<Double> Ys = new ArrayList<Double>();
				for (int i = 0; i < RT1.size(); i++) {
					double X = RT1.getValue("X", i);
					double Y = RT1.getValue("Y", i);
					Xs.add(X);
					Ys.add(Y);
				}
				matrix.add(Xs);
				matrix.add(Ys);
				Values_RT.add(matrix);
			}
			IJ.selectWindow("Results"); 
			IJ.run("Close");
			return Values_RT;
	}
	
	public RoiManager track_wells(ImagePlus input) {
		// Extract only one image as the wells do not move across time 
		ImagePlus imp_ch = input.crop("whole-slice");
		imp_ch.setTitle("Segmented channels");
		IJ.run(imp_ch, "Find Edges", "stack");
		IJ.setAutoThreshold(imp_ch, "Default dark no-reset");
		Prefs.blackBackground = false;
		IJ.run(imp_ch, "Convert to Mask", "method=Default background=Dark calculate");
		IJ.run(imp_ch, "Fill Holes", "stack");
		// Close the Threshold window and its output
		//IJ.run("Close", "Threshold");
		//IJ.resetThreshold(input);
		
		// Analyse the particles TODO need fine tuning of parameters
		IJ.run(imp_ch, "Analyze Particles...", "size=1000-100000 circularity=0.00-0.50 show=[Overlay Masks] exclude clear overlay add");
		//imp_ch.show();
		
		// Work with rois
		IJ.run("From ROI Manager", "");
		RoiManager rm = RoiManager.getInstance();
		return rm;
	}
	
	void draw_wells(RoiManager rm, ArrayList<ArrayList<ArrayList<Double>>> cells_position) {
		int n = rm.getCount();
		rm.runCommand("List"); 
		ResultsTable RT1 = ResultsTable.getResultsTable("Overlay Elements of Result of p1_stabilised_cells.tif"); //TODO: Change the title
		double wells_width = RT1.getValue("Width", 0);	// Retrieve width of wells
		double wells_height = RT1.getValue("Height", 0); // Retrieve height of wells
		ArrayList<Double> Xs_wells = new ArrayList<Double>();
		ArrayList<Double> Ys_wells = new ArrayList<Double>();
		for (int i = 0; i < RT1.size(); i++) {
			double X = RT1.getValue("X", i); // Retrieve the position in the upper-left corner of the rectangle around the wells
			double Y = RT1.getValue("Y", i);
			Xs_wells.add(X);
			Ys_wells.add(Y);
		}
		ArrayList<ArrayList<Double>> first_frame = cells_position.get(0);
		ArrayList<Double> Xs_cells = first_frame.get(0);
		ArrayList<Double> Ys_cells = first_frame.get(1);
		for(int i=0;i<=n-1;i++) {
			rm.select(i);
			boolean isInside = false;
			// Check if any cells are inside any wells
			for (int k=0; k<=Xs_cells.size()-1; k++) {
				if (Xs_cells.get(k) > Xs_wells.get(i) && Xs_cells.get(k) < Xs_wells.get(i) + wells_width & Ys_cells.get(k) > Ys_wells.get(i) & Ys_cells.get(k) < Ys_wells.get(i) + wells_height) {
					isInside = true;
				}
			}
			if (isInside) {
				rm.runCommand("Set Color", "green"); // If inside -> color = green
			}
			else {
				rm.runCommand("Set Color", "red"); // if not inside -> color = red
			}
		}
		rm.runCommand("Show all");
	}
	
	public ImagePlus remove_background(RoiManager rm, ImagePlus wells, ImagePlus cells) {
		Roi[] listRois  = rm.getRoisAsArray();
		Overlay overlay = new Overlay(); 
		for(int i=0;i<=listRois.length-1;i++) {
			// Gather the overlays from the roi manager and shrink them  
			overlay.add(RoiEnlarger.enlarge(listRois[i], -3)); 
		}
		wells.setOverlay(overlay);
		// Extract the mask from the overlays
		ByteProcessor mask = wells.createRoiMask();
		//input.show(); 
		// Convert the mask to a binary image
		ImageStack stack = new ImageStack(); 
		stack.addSlice(mask); 
		ImagePlus impMask = new ImagePlus("Masks", stack); 
		// Convert the image to 16 bits, and stretch to maximal values
		IJ.run(impMask, "16-bit", "");
		IJ.run(impMask, "Multiply...", "value=258");
		// Invert such that the channels are black (value=0) and surroundings are white
		IJ.run(impMask, "Invert", "");
		//impMask.show();
		// Remove background to original image
		ImagePlus imp3 = ImageCalculator.run(cells, impMask, "Subtract create stack");
		imp3.show(); 
		return imp3; 
	}

	private void draw(Overlay overlay, ArrayList<Spot> spots[]) {
		int nt = spots.length;
		for (int t = 0; t < nt; t++)
			for (Spot spot : spots[t])
				spot.draw(overlay);
	}

	// Array of ArrayList of Spot, first array is per frame, the second is the number of spots per array? 
	private ArrayList<Spot>[] filter(ImagePlus dog, ArrayList<Spot> spots[], double threshold) {
		int nt = spots.length;
		// Create an array of size # frames, containing ArrayList of Spot 
		ArrayList<Spot> out[] = new ArrayList[nt];
		for (int t = 0; t < nt; t++) {
			// Initialise each row of the array by an empty ArrayList of Spot 
			out[t] = new ArrayList<Spot>();
			for (Spot spot : spots[t]) {
				// Add the spot "overlay" to the image? 
				dog.setPosition(1, 1, t + 1);
				double value = dog.getProcessor().getPixelValue(spot.x, spot.y);
				// If the spot from the list is higher than the threshold, add it 
				if (value > threshold)
					out[t].add(spot);
			}
		}
		return out;
	}

	private ImagePlus dog(ImagePlus imp, double sigma) {
		ImagePlus g1 = imp.duplicate();
		ImagePlus g2 = imp.duplicate();
		IJ.run(g1, "Gaussian Blur...", "sigma=" + sigma + " stack");
		IJ.run(g2, "Gaussian Blur...", "sigma=" + (Math.sqrt(2) * sigma) + " stack");
		ImagePlus dog = ImageCalculator.run(g1, g2, "Subtract create stack");
		dog.setTitle("DOG of input image");
		dog.show();
		return dog;
	}


	private ArrayList<Spot>[] localMax(ImagePlus imp) {
		int nt = imp.getNFrames();
		int nx = imp.getWidth();
		int ny = imp.getHeight();
		// Create an array of Spots with size = nb. of frames 
		// Figure out what a Spots is
		ArrayList<Spot> spots[] = new ArrayList[nt];
		for (int t = 0; t < nt; t++) {
			imp.setPosition(1, 1, t + 1);
			ImageProcessor ip = imp.getProcessor();
			spots[t] = new ArrayList<Spot>();
			// Loop over the image (except the border apparently) 
			for (int x = 1; x < nx - 1; x++) {
				for (int y = 1; y < ny - 1; y++) {
					double v = ip.getPixelValue(x, y);
					double max = -1;
					// Loop on a 3x3 square around the pixel (x,y), and get the max of the neighbouring pixels
					for (int k = -1; k <= 1; k++)
						for (int l = -1; l <= 1; l++)
							max = Math.max(max, ip.getPixelValue(x + k, y + l));
					// If (x,y) pixel is the brightest is a neighbourhood of 3x3 pixels, add it as a spot 
					if (v == max)
						spots[t].add(new Spot(x, y, t, v));
				}
			}
		}
		return spots;
	}

}