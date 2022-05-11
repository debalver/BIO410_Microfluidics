import java.util.ArrayList;
import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import ij.gui.Overlay;
import ij.plugin.ImageCalculator;
import ij.plugin.PlugIn;
import ij.process.ImageProcessor;

public class Draw_Tracking implements PlugIn {

	public void run(String arg) {
		double sigma = 3;
		double threshold = 10;
		double lambda = 0.5;
		double threshold_2 = 1e-2;
		boolean isNearest = false;
		ImagePlus imp1 = IJ.getImage();
		ImagePlus imp = imp1.duplicate();
		imp1.close();
		IJ.run(imp, "Size...", "width=600 height=400 depth=260 average interpolation=Bilinear");
		ImagePlus imp_final = imp.duplicate();
		imp.show();
		IJ.setThreshold(235, 255);
		IJ.run(imp, "Make Binary", "method=Default background=Default black");
		imp.setTitle("binary");
		IJ.selectWindow("binary");
		IJ.run("32-bit", "");
		int nt = imp.getNFrames();
		ImagePlus dog = dog(imp, sigma);
		ArrayList<Spot> localmax[] = localMax(dog);
		ArrayList<Spot> spots[] = filter(dog, localmax, threshold);
		double width = imp.getWidth();
		double height = imp.getHeight();
		double d_max = Math.sqrt(width * width + height * height);
		double f_max = getFmax(spots, nt, dog);
		for (int t = 0; t < nt - 1; t++) {
			for (Spot current : spots[t]) {
				for (Spot next : spots[t+1]) {
					if (isNearest) {
						current.store_neighbors(next);
					}
					else {
						if (current.distance(next, lambda, d_max, f_max) < threshold_2) {
							current.link(next);
						}
					}
				}
				if (isNearest) {
					current.nearest_neighbors(lambda, d_max, f_max, threshold_2);
				}
			}
		}
		Overlay overlay = new Overlay();
		draw(overlay, spots);
		imp_final.setOverlay(overlay);
	}

	private void draw(Overlay overlay, ArrayList<Spot> spots[]) {
		int nt = spots.length;
		for (int t = 0; t < nt; t++)
			for (Spot spot : spots[t])
				spot.draw(overlay);
	}

	private ArrayList<Spot>[] filter(ImagePlus dog, ArrayList<Spot> spots[], double threshold) {
		int nt = spots.length;
		ArrayList<Spot> out[] = new Spots[nt];
		for (int t = 0; t < nt; t++) {
			out[t] = new Spots();
			for (Spot spot : spots[t]) {
				dog.setPosition(1, 1, t + 1);
				double value = dog.getProcessor().getPixelValue(spot.x, spot.y);
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
		dog.show();
		return dog;
	}


	private Spots[] localMax(ImagePlus imp) {
		int nt = imp.getNFrames();
		int nx = imp.getWidth();
		int ny = imp.getHeight();
		Spots spots[] = new Spots[nt];
		for (int t = 0; t < nt; t++) {
			imp.setPosition(1, 1, t + 1);
			ImageProcessor ip = imp.getProcessor();
			spots[t] = new Spots();
			for (int x = 1; x < nx - 1; x++) {
				for (int y = 1; y < ny - 1; y++) {
					double v = ip.getPixelValue(x, y);
					double max = -1;
					for (int k = -1; k <= 1; k++)
						for (int l = -1; l <= 1; l++)
							max = Math.max(max, ip.getPixelValue(x + k, y + l));
					if (v == max)
						spots[t].add(new Spot(x, y, t, v));
				}
			}
		}
		return spots;
	}

	
	private double getFmax(ArrayList<Spot> spots[], int nt, ImagePlus dog) {
		double F_max = 0;
		for (int t = 0; t < nt; t++) {
			for (Spot current : spots[t]) {
				double temp = dog.getProcessor().getPixelValue(current.x, current.y);
				if (F_max < temp) {
					F_max = temp;
				}
			}
		}
		return F_max;
	}
}