import java.awt.Color;
import java.util.ArrayList;

import ij.gui.Line;
import ij.gui.OvalRoi;
import ij.gui.Overlay;
import ij.gui.Roi;

public class Cell {
	public int x;
	public int y;
	public int t;
	public double intensity; 
	private Cell next = null;
	private Color color;
	private ArrayList<Cell> list_neigh = new ArrayList<Cell>(); 

	public Cell(int x, int y, int t, double i) {
		this.x = x;
		this.y = y;
		this.t = t;
		this.intensity = i; 
		color = Color.getHSBColor((float)Math.random(), 1f, 1f);
		color = new Color(color.getRed(), color.getGreen(), color.getBlue(), 120);
	}

	public double distance(Cell spot, double l, double d_max, double s_max) {
		// Basic Euclidean distance  
		double dx = x - spot.x;
		double dy = y - spot.y;
		double dist =  Math.sqrt(dx * dx + dy * dy);
		// Get the intensity difference 
		double similarity = Math.abs(intensity-spot.intensity);
		// More complex cost function 
		return (1-l)*dist/d_max + l*similarity/s_max; 
	}

	public void draw(Overlay overlay) {
		double xp = x + 0.5;
		double yp = y + 0.5;
		int radius = 5;
		OvalRoi roi = new OvalRoi(xp - radius, yp - radius, 2 * radius, 2 * radius);
		roi.setPosition(t+1); // display roi in one frqme
		roi.setStrokeColor(new Color(255, 0, 0, 120));
		roi.setStrokeWidth(1);
		overlay.add(roi);
		if (next != null) {
			Line line = new Line(x, y, next.x, next.y);
			line.setStrokeColor(color);
			line.setStrokeWidth(2);
			overlay.add(line);
		}
		
		//TextRoi text = new TextRoi(x, y-10, "" + value);
		//text.setPosition(t+1);
		//overlay.add(text);
	}
	
	public void link(Cell a) {
		if (a == null)
			return;
		a.next = this;
		a.color = this.color;
	}
	
	public void add_neighbour(Cell a) {
		if (a != null) 
			this.list_neigh.add(a); 
	}
	
	public void get_nearest_neighbour(double lambda, double d_max, double s_max, double thres) {
		double c_min = Double.POSITIVE_INFINITY;
		Cell best = null; 
		// Iterate on the neighbours
		for(Cell neigh : this.list_neigh) {
			double c = this.distance(neigh, lambda, d_max, s_max);
			 // If neighbour cost is smaller than the previous ones and below the threshold, link it 
			if (c < c_min && c < thres)  
				best = neigh;
				c_min = c; 
		}
		this.link(best); 
	}
	
	public String toString() {
		return "(" + x + ", " + y + ", " + t + ")";
	}
}