import java.util.ArrayList;
import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import ij.gui.Overlay;
import ij.plugin.ImageCalculator;
import ij.plugin.PlugIn;
import ij.process.ImageProcessor;
import ij.ImagePlus;
import ij.plugin.PlugIn;
import ij.*;
import ij.plugin.filter.*;
import ij.plugin.frame.RoiManager;
import ij.process.*;
import ij.gui.*;
import ij.measure.*;


public class Draw_Channels implements PlugIn {
	public void run(String arg) {
		ImagePlus imp = IJ.getImage();
		IJ.run("Make Substack...", "slices=1");
		ImagePlus sub = IJ.getImage();
		IJ.run("Convert to Mask");
		IJ.run("Analyze Particles...", "size=1000-100000 show=Overlay exclude clear overlay add");
		IJ.run("From ROI Manager", "");
		RoiManager rm = RoiManager.getInstance();
		int n = rm.getCount();
		for(int i=0;i<=n-1;i++) {
			rm.select(i);
			if (i == 5) { //C'est juste pour voir que ça fonctionne mais il faudra avoir une autre condition ici
				rm.runCommand("Set Color", "red");
			}
			else {
				rm.runCommand("Set Color", "green");
			}
		}
		sub.close();
		rm.runCommand(imp, "Show All without labels");
	}
}
	