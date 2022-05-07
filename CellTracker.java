import java.util.ArrayList;
import ij.IJ;
import ij.ImagePlus;
import ij.gui.Overlay;
import ij.plugin.ImageCalculator;
import ij.plugin.PlugIn;
import ij.process.ImageProcessor;
// No need to import Spot if it is in the same folder

public class CellTracker implements PlugIn {
	@Override
	public void run(String arg) {
		
		ImagePlus imp = IJ.getImage();
		ImagePlus imp_ch = imp.crop("slice=2"); 
		//IJ.run("Duplicate...", "duplicate slices=2");
		IJ.run(imp_ch, "Find Edges", "stack");
		//IJ.run("Threshold...");
		IJ.setAutoThreshold(imp_ch, "Minimum dark");
		//IJ.run("Threshold...");
		IJ.run(imp_ch, "Analyze Particles...", "size=150-Infinity circularity=0.00-0.13 show=[Overlay Masks] clear add stack");
		imp_ch.show();
		
		
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