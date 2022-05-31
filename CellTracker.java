import java.util.ArrayList;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.GenericDialog;
import ij.gui.Overlay;
import ij.gui.Roi;
import ij.plugin.ImageCalculator;
import ij.plugin.PlugIn;
import ij.plugin.RoiEnlarger;
import ij.plugin.frame.RoiManager;
import ij.process.ByteProcessor;
import ij.process.ImageStatistics;
import ij.Prefs; 
import ij.WindowManager; 
import ij.gui.Overlay; 
import ij.measure.ResultsTable;
import javax.swing.*; 

// By Quentin Devaud & Laurent Gürtler - May 2022 

public class CellTracker implements PlugIn {
	@Override
	public void run(String arg) {
		
		////////////////////////////////////////// Part for dialog box /////////////////////////////////////////////////////
		
		ImagePlus original = null;
		ImagePlus aligned_stack_wells = null; 
		ImagePlus aligned_stack_cells = null;
		double slice_wells = 0; 
		double slice_cells = 0; 
		boolean align = false; 
		
		GenericDialog gui = new GenericDialog("Stabilised"); 
		String[] choices = {"Yes", "No"};
		gui.addChoice("Do you have stabilised images?", choices, "No"); 
		gui.showDialog();
		String choice = gui.getNextChoice();
		if (gui.wasOKed()) {
			if (choice=="No") {
				GenericDialog gui2 = new GenericDialog("Choose file");
				gui2.addNumericField("Slice with the sharpest well: ", 2); 
				gui2.addNumericField("Slice with cells having the best contrast: ", 1);
				gui2.addFileField("Choose the images without stabilisation: ", "/path/image");
				gui2.showDialog();
				align = true; 
				slice_wells = gui2.getNextNumber();
				slice_cells = gui2.getNextNumber(); 
				String path = gui2.getNextString(); 
				original = IJ.openImage(path);
				
			} else {
				GenericDialog gui3 = new GenericDialog("Choose files");
				gui3.addFileField("Choose the images with the sharpest wells: ", "/path/image_wells");
				gui3.addFileField("Choose the images with cells having the best contrast: ", "/path/image_cells");
				gui3.showDialog();
				
				if (gui3.wasOKed()) {
					align = false; 
					String path_wells = gui3.getNextString();
					String path_cells = gui3.getNextString(); 
					aligned_stack_wells = IJ.openImage(path_wells);
					aligned_stack_cells = IJ.openImage(path_cells);
				}
			}
		}
		

		////////////////////////////////////////// Part for aligning the images /////////////////////////////////////////////////////
		
		// During development we directly work on the aligned images
		if (align) {
			// Meaning of the numbers in Duplicator:
			// firstC, lastC, firstZ, lastZ, firstT, lastT
			// For experiment 1, the slice 2 is better for the channel detection
			int bestZ_wells = (int)slice_wells;
			int bestZ_cells = (int)slice_cells; 
			int lastT = original.getNFrames();
			// Check if the slices for cells and wells are the same, if yes: does less computations
			if (bestZ_wells != bestZ_cells) {
				// Extract only the  most useful slices
				aligned_stack_wells = new ij.plugin.Duplicator().run(original, 1, 1, bestZ_wells, bestZ_wells, 1, lastT); 
				aligned_stack_cells = new ij.plugin.Duplicator().run(original, 1, 1, bestZ_cells, bestZ_cells, 1, lastT);
				// Save memory and "delete" the original
				original = null; 
				// Align the images across time 
				message("Starting the alignment of wells..."); 
				IJ.run(aligned_stack_wells, "StackReg ", "transformation=Affine");
				message("Alignment of wells done, moving to cells alignment..."); 
				IJ.run(aligned_stack_cells, "StackReg ", "transformation=Affine");
				message("Alignment done, moving to tracking!"); 
			} else {
				aligned_stack_wells = new ij.plugin.Duplicator().run(original, 1, 1, bestZ_wells, bestZ_wells, 1, lastT); 
				// Save memory and "delete" the original
				original = null; 
				// Align the images across time 
				message("Aligning the images of the wells and cells...");
				IJ.run(aligned_stack_wells, "StackReg ", "transformation=Affine");
				aligned_stack_cells = aligned_stack_wells; 
			}
		} 
		 
		aligned_stack_wells.setTitle("cells_microfluidics"); 
		aligned_stack_cells.setTitle("cells_microfluidics");
		
		
		////////////////////////////////////////// Parts for tracking the cells and wells /////////////////////////////////////////////////////
		
		RoiManager rm = track_wells(aligned_stack_wells);
		// Returns the masked image as well as the shrank ROIs 
		Object[] temp  = remove_background(rm, aligned_stack_wells.duplicate(), aligned_stack_cells.duplicate());
		ImagePlus cells_only = (ImagePlus)temp[0]; 
		RoiManager s_rm = (RoiManager)temp[1]; 
		Object[] stats = mean_Roi_intensity(s_rm, cells_only);
		ArrayList<ArrayList<ArrayList<Double>>> cells_position = track_cells(cells_only, (double)stats[0], (double)stats[1]); // Position of the cells in terms of numbers of frames, X and Y, number of cells per frame
		draw_wells(s_rm, cells_position); // Draw color if cell is inside well 
		ArrayList<Cell> cells[] = create_cells(cells_position, aligned_stack_cells); 
		
		////////////////////////////////////////// Parts for linking the cells /////////////////////////////////////////////////////
		
		// Some parameters for computing the cost function for linking
		int nt = aligned_stack_cells.getNFrames();
		int nx = aligned_stack_cells.getWidth();
		int ny = aligned_stack_cells.getHeight();
		
		double d_max = Math.sqrt(ny*ny+nx*nx); 
		double s_max = Math.pow(2, aligned_stack_cells.getBitDepth()); 
		double lambda = 0.2; 
		double thres = 5.5e-2; 
		boolean nearest = true; 
		
		Overlay overlay = new Overlay();
		// Draw the trajectories for the previous x points
		int x = 10; 
		
		// Loop on the frames: from the second to the last 
		for (int t = 1; t < nt; t++) {
			// Loop on the spots of the frame t
			for (Cell current : cells[t]) {
				// Loop on the spots of the next frame 
				for (Cell previous : cells[t-1]) {
					// Nearest neighbour linking: first store all the potential neighbours
					if(nearest)
						current.add_neighbour(previous);
					// Simple linking: link only if the cost is below a threshold 
					else
						if (current.distance(previous, lambda, d_max, s_max) < thres)
							current.link_previous(previous);
				// Then get the nearest neighbour among all the potential neighbours 
				}
				if(nearest) {
					current.get_nearest_neighbour(lambda, d_max, s_max, thres);}
				
				// Draw the trajectory line for the last x points
				current.draw_line(overlay, x, t); 
			}
		}
		
		// Draw the red square for the cell locations 
		draw_location(overlay, cells);
		// Display the trajectories and red squares
		aligned_stack_cells.setOverlay(overlay);
		aligned_stack_cells.show(); 
		// Display the wells above the cells 
		s_rm.runCommand("Show all");
	
	}
	
	////////////////////////////////////////// All useful functions /////////////////////////////////////////////////////
		

	
	/**
	 * A function that tracks the cells present in an image. 
	 * 
	 * @param imp: The image used to track the cells 
	 * @return Values_RT: A triple array containing first the frames, then X, then Y positions of the cells
	 */
	public ArrayList<ArrayList<ArrayList<Double>>> track_cells(ImagePlus imp, double max, double std) {
			
			// Get a mask with only the cells based on standard deviation and max value of the well
			imp.setDisplayRange(max-std*1.2, max+100);
			IJ.run(imp, "Apply LUT", "stack");
			IJ.setAutoThreshold(imp, "Default no-reset");
			Prefs.blackBackground = false;
			IJ.run(imp, "Convert to Mask", "method=Default background=Light calculate");
			// Blur the image to remove non relevant particles
			IJ.run(imp, "Gaussian Blur...", "sigma=6 stack");
			int nt = imp.getNFrames();
			// Triple array containing the cells positions. First dim. Frame, Second dim. X , Thrid dim. Y
			ArrayList<ArrayList<ArrayList<Double>>> Values_RT = new ArrayList<ArrayList<ArrayList<Double>>>();
			for (int t = 0; t < nt; t++) { // loop over frames
				ArrayList<ArrayList<Double>> matrix = new ArrayList<ArrayList<Double>>();
				imp.setSlice(t); // select slice
				// Find the maxima 
				IJ.run(imp, "Find Maxima...", "prominence=10 output=List");
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
	
	
	/**
	 * A function which tracks the wells present in the image. 
	 * 
	 * @param input: The image analysed for the finding the wells. 
	 * @return rm: A RoiManager with the shape and coordinates of the wells. 
	 */
	public RoiManager track_wells(ImagePlus input) {
		// Extract only one image as the wells do not move across time 
		ImagePlus imp_ch = input.crop("whole-slice");
		imp_ch.setTitle("Segmented channels");
		// Start by finding the edges and set a threshold and convert to mask
		IJ.run(imp_ch, "Find Edges", "stack");
		IJ.setAutoThreshold(imp_ch, "Default dark no-reset");
		Prefs.blackBackground = false;
		IJ.run(imp_ch, "Convert to Mask", "method=Default background=Dark calculate");
		// Fill the wells 
		IJ.run(imp_ch, "Fill Holes", "stack");
		// Analyse only the large particles which are not too round
		IJ.run(imp_ch, "Analyze Particles...", "size=1000-100000 circularity=0.00-0.50 show=[Overlay Masks] exclude clear overlay add");
		// Get the ROI manager and its results table - containing width, and height of the ROIs 
		RoiManager rm = RoiManager.getInstance();
		rm.runCommand(imp_ch,"List");
		
		return rm;
	}
	
	/**
	 * This function modifies the colour of the well's ROI to green if a cell is inside. 
	 * 
	 * @param rm: A RoiManager containing the wells 
	 * @param cells_position: A triple array containing the frame, X, and Y positions of the cells 
	 */
	void draw_wells(RoiManager rm, ArrayList<ArrayList<ArrayList<Double>>> cells_position) {
		int n = rm.getCount();
		ResultsTable RT_wells = ResultsTable.getResultsTable("Overlay Elements of Segmented channels");
		double wells_width = RT_wells.getValue("Width", 0);	// Retrieve width of wells
		double wells_height = RT_wells.getValue("Height", 0); // Retrieve height of wells
		ArrayList<Double> Xs_wells = new ArrayList<Double>();
		ArrayList<Double> Ys_wells = new ArrayList<Double>();
		for (int i = 0; i < RT_wells.size(); i++) {
			double X = RT_wells.getValue("X", i); // Retrieve the position in the upper-left corner of the rectangle around the wells
			double Y = RT_wells.getValue("Y", i);
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
		IJ.selectWindow("Overlay Elements of Segmented channels");
		IJ.run("Close");
		
	}
	
	
	/**
	 * This function transforms a set of coordinates to an object of class Cell and returns an array of Cell.  
	 * 
	 * @param cells_position: A triple array containing the frame, X, and Y positions of the cells
	 * @param aligned_cells: The image of the aligned cells 
	 * @return cells: An array list of array of Cell containing the cells of each frame. 
	 */
	public ArrayList<Cell>[] create_cells(ArrayList<ArrayList<ArrayList<Double>>> cells_position, ImagePlus aligned_cells) {
		// Get the number of frames in the image
		int nb_f = cells_position.size();
		ArrayList<Cell> cells[] = new ArrayList[nb_f];
		// Iterate on the frames 
		for(int t=0; t<=nb_f-1; t++) {
			// Extract the X and Y positions for the frame t  
			ArrayList<Double> Xs_cells = cells_position.get(t).get(0);
			ArrayList<Double> Ys_cells = cells_position.get(t).get(1);
			cells[t] = new ArrayList<Cell>();
			// Iterate over the X and Y positions
			for (int k=0; k<=Xs_cells.size()-1; k++) {
				// Convert back the coordinate to int 
				int x = Xs_cells.get(k).intValue(); 
				int y = Ys_cells.get(k).intValue(); 
				// Get the pixel average pixel value around the cell 
				double i = get_avg_pixel_value(x, y, aligned_cells); 
				// Create the spot and add it to the list
				cells[t].add(new Cell(x, y, t, i));
			}
		}
		
		return cells; 
	}
	
	/**
	 * This function computes the average pixel value around the pixel x,y in the image imp. 
	 * 
	 * @param x
	 * @param y
	 * @param imp
	 * @return
	 */
	public double get_avg_pixel_value(int x, int y, ImagePlus imp) {
		double intensity = 0;
		// Iterate on a 3x3 square around the pixel and get the pixel values 
		for (int k = -1; k <= 1; k++)
			for (int l = -1; l <= 1; l++)
				intensity += new Double(imp.getPixel(x + k, y + l)[0]); // [0] should get the gray scale value 
		// And average the values inside of this square	
		return intensity/9;
	}
	
	/**
	 * This function sets to 0 all pixels value of the image "cells", if the pixels are located outside of the wells in the ROI manager rm.  
	 * 
	 * @param rm: The ROI manager containing the shape and locations of the wells. 
	 * @param wells: The image with the sharpest wells 
	 * @param cells: The image with the cells having the best contrast with respect to background
	 * @return : Returns the masked image as well as the shrunk ROIs (which better fit the wells) 
	 */
	public Object[] remove_background(RoiManager rm, ImagePlus wells, ImagePlus cells) {
		wells.setTitle("cells_microfluicis"); 
		cells.setTitle("cells_microfluicis"); 
		
		Roi[] listRois  = rm.getRoisAsArray();
		// New ROI manager which stores the shrunk ROIs 
		RoiManager s_rm = new RoiManager(true); 
		Overlay overlay = new Overlay(); 
		for(int i=0;i<=listRois.length-1;i++) {
			// Gather the overlays from the roi manager and shrink them
			Roi small_roi = RoiEnlarger.enlarge(listRois[i], -3);
			overlay.add(small_roi);
			s_rm.addRoi(small_roi); 
		}
		wells.setOverlay(overlay);
		// Extract the mask from the overlays
		ByteProcessor mask = wells.createRoiMask();
		// Convert the mask to a binary image
		ImageStack stack = new ImageStack(); 
		stack.addSlice(mask); 
		ImagePlus impMask = new ImagePlus("Masks", stack); 
		// Convert the image to 16 bits, and stretch to maximal values
		IJ.run(impMask, "16-bit", "");
		IJ.run(impMask, "Multiply...", "value=258");
		// Invert such that the channels are black (value=0) and surroundings are white
		IJ.run(impMask, "Invert", "");
		// Remove background to original image
		ImagePlus imp3 = ImageCalculator.run(cells, impMask, "Subtract create stack");

		return new Object[] {imp3, s_rm}; 
	}

	private void draw_location(Overlay overlay, ArrayList<Cell> cells[]) {
		int nt = cells.length;
		for (int t = 0; t < nt; t++)
			for (Cell cell : cells[t])
				cell.draw(overlay);
	}
	
	public static void message(String msg) {
        JOptionPane.showMessageDialog(null,
               msg,
               "PopUp Dialog",
               JOptionPane.INFORMATION_MESSAGE);
   }
	
	/**
	 * Computes the max value and standard deviation of the brightest well.    
	 * 
	 * @param rm: The ROI manager containing the intensity inside the wells.
	 * @param imp: The image on which we compute the statistics 
	 * @return : Returns the mean intensity 
	 */
	public Object[] mean_Roi_intensity(RoiManager rm, ImagePlus imp) {
		double max = Double.NEGATIVE_INFINITY; 
		double std = Double.POSITIVE_INFINITY; 
		Roi[] rois = rm.getRoisAsArray(); 
		// Iterate on ROIs 
		for(int i=0; i<rois.length ;i++) {
			// Add the ROI of interest and get its statistics 
			imp.setRoi(rois[i]);
			ImageStatistics stats = imp.getStatistics(); 
			// Extract the highest value and standard deviation
			if (stats.max > max) {
				max = stats.max; 
				std = stats.stdDev; 
			}
				
			imp.resetRoi(); 
		}
		return new Object[] {max, std}; 
	}
		
	// End class
}
