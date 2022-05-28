import java.awt.Color;
import java.util.ArrayList;

import ij.gui.Line;
import ij.gui.OvalRoi;
import ij.gui.Overlay;
import ij.gui.PointRoi;
import ij.gui.PolygonRoi;
import ij.gui.Roi;

public class Cell {
	public int x;
	public int y;
	public int t;
	public double intensity; 
	private Cell next = null;
	public Cell previous = null; 
	public Color color;
	private ArrayList<Cell> list_neigh = new ArrayList<Cell>(); 

	public Cell(int x, int y, int t, double i) {
		this.x = x;
		this.y = y;
		this.t = t;
		this.intensity = i; 
		color = Color.getHSBColor((float)Math.random(), 1f, 1f);
		color = new Color(color.getRed(), color.getGreen(), color.getBlue(), 120);
	}

	public double distance(Cell cell, double l, double d_max, double s_max) {
		// Basic Euclidean distance  
		double dx = x - cell.x;
		double dy = y - cell.y;
		double dist =  Math.sqrt(dx * dx + dy * dy);
		// Get the intensity difference 
		double similarity = Math.abs(intensity-cell.intensity);
		// More complex cost function 
		return (1-l)*dist/d_max + l*similarity/s_max; 
	}

	public void draw(Overlay overlay) {
		double xp = x;
		double yp = y;
		int radius = 8;
		// Draw the ROI locating the cell (red circle) 
		//OvalRoi roi = new OvalRoi(xp - radius, yp - radius, 2 * radius, 2 * radius);
		// Draw the ROI locating the cell (red square) 
		Roi roi = new Roi(xp, yp, radius, radius); 
		// Only display the ROI in current frame 
		roi.setPosition(t+1); // display roi in one frqme
		roi.setStrokeColor(new Color(235, 0, 0, 150));
		roi.setStrokeWidth(1);
		overlay.add(roi);
	}
	
	public void draw_line(Overlay overlay, int x, int t) {
		
		Cell curr = this; 
		Cell prev = this.previous;
		Color color = this.color; 
		
		for(int p = 0; p < x; p++) {
			if(curr != null & prev != null) {
				Line line = new Line(prev.x, prev.y, curr.x, curr.y);
				line.setStrokeColor(color);
				line.setStrokeWidth(2);
				// Only draw the line for the given frame 
				line.setPosition(t+1);
				overlay.add(line);
				// Move back in time by one frame 
				curr = prev; 
				prev = prev.previous; 
			}
		}
		
	}
	
	public void link_previous(Cell a) {
		// Link the current cell to the previous one 
		if (a == null)
			return;
		//a.next = this;
		//a.color = this.color;
		this.previous = a; 
		this.color = a.color; 
	}
	
	public void link_next(Cell a) {
		// Link the current cell to the next one 
		if (a == null)
			return;
		// Check if the current cell already has a next cell
		// If yes, a division took place, then change the next cell colour 
		if(this.next != null)
			a.color = Color.getHSBColor((float)Math.random(), 1f, 1f);
		this.next = a;  
	}
	
	public void add_neighbour(Cell a) {
		if (a != null) 
			this.list_neigh.add(a); 
	}
	
	/**
	 * A function that finds the nearest previous neighbour of the given cell.
	 * 
	 * @param lambda: Adjusts the cost function towards more weight for the distance or intensity 
	 * @param d_max: The maximum distance available on the image 
	 * @param s_max: The maximum intensity available on the image 
	 * @param thres: The minimum threshold (in terms of cost) to link two cells 
	 */
	public void get_nearest_neighbour(double lambda, double d_max, double s_max, double thres) {
		double c_min = Double.POSITIVE_INFINITY;
		Cell best = null; 
		// Iterate on the neighbours
		for(Cell neigh : this.list_neigh) {
			double c = this.distance(neigh, lambda, d_max, s_max);
			 // If neighbour cost is smaller than the previous ones and below the threshold, link it 
			if ( ((c < c_min)==true) && ((c < thres)==true) ) {
				best = neigh;
				c_min = c; 
			}
		}
		this.link_previous(best);
		if(best != null)
			best.link_next(this); 
	}
	
	public String toString() {
		return "(" + x + ", " + y + ", " + t + ")";
	}
}