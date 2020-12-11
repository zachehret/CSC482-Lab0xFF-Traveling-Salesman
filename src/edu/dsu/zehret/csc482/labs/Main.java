package edu.dsu.zehret.csc482.labs;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Random;

public class Main {

	private static int VERTEX_COUNT = 10;
	private static Vertices vertices = new Vertices(VERTEX_COUNT);
	public static void main(String[] args) {
		Random ran = new Random();

		int i = 0;
		while(i < 20) {
			VERTEX_COUNT = 5 + 1*i; //10...20...30...40...50......
			vertices = new Vertices(VERTEX_COUNT);
			System.out.println("#################################");
			System.out.println(VERTEX_COUNT + " vertices");
			int[][] m = generateRandomCostMatrix(VERTEX_COUNT, 50, ran);
			long startTime = System.currentTimeMillis();
			int greedy = vertices.greedyAlgorithmCost(m);
			long endTime = System.currentTimeMillis();
			System.out.println("Greedy: " + greedy + " (" + ((endTime - startTime == 0) ? "<1" : (endTime - startTime)) + "ms)");
			
			startTime = System.currentTimeMillis();
			int bruteForce = vertices.bruteForceAlgorithmCost(m);
			endTime = System.currentTimeMillis();
			System.out.println("Brute Force: " + bruteForce + " (" + ((endTime - startTime == 0) ? "<1" : (endTime - startTime)) + "ms)");
			
			i++;
		}
	}
	
	/**
	 *
	 * @param n Number of vertices
	 * @param c Max Cost
	 * @return Matrix
	 * 
	 */
	private static int[][] generateRandomCostMatrix(int n, int c, Random ran) {
		int[][] m = new int[n][n];
		for(int i = 0; i < n; i++) {
			int j = i;
			while(j > 0) {
				//System.out.println(i + "+" + j + "=" + (i+j));
				//System.out.println(i + "-" + j + "=" + (i-j));
				int r = ran.nextInt(c) + 1;
				if(i+j < n) {
					m[i+j][i+j] = r;
					m[i+j][i] = r;
					m[i][i+j] = r;
				}
				if(i-j >= 0) {
					m[i-j][i-j] = r;
					m[i-j][i] = r;
					m[i][i-j] = r;
				}				
				j--;
			}
		}
		for(int i = 0; i < n; i++) {
			m[i][i] = 0;
		}
		return m;
	}
	
	/**
	 * 
	 * @param n Number of vertices
	 * @param maxX max x coordinate
	 * @param maxY max y coordinate
	 * @param ran random
	 * @return Cost Matrix
	 */
	private static int[][] generateRandomEuclideanCostMatrix(int n, int maxX, int maxY, Random ran) {
		Coordinate[] coordinates = new Coordinate[n];
		for(int i = 0; i < coordinates.length; i++) {
			//Generate random coordinates
			coordinates[i] = new Coordinate(maxX, maxY, ran);
		}
		int[][] m = new int[n][n];
		for(int i = 0; i < coordinates.length; i++) {
			for(int j = 0; j < coordinates.length; j++) {
				m[i][j] = (int) coordinates[i].distanceFrom(coordinates[j]);
			}
		}
		return m;
	}
	
	/**
	 * 
	 * @param n Number of vertices
	 * @param r Radius
	 * @return Cost Matrix
	 */
	private static float[][] generateRandomCircularGraphCostMatrix(int n, int r) {
		//Distance per vertex = circumference / n
		float distancePerVertice = (float) ((2f * Math.PI * r) / n) ;
		float[][] m = new float[n][n];
		
		for(int i = 0; i < n; i++) {
			for(int j = 0; j < n; j++) {
				if(i == j)
					m[i][j] = 0;
				else
					m[i][j] = (float) distancePerVertice;
			}
		}
		return m;
	}
}
class Coordinate {
	public float x, y;
	public Coordinate(float maxX, float maxY, Random ran) {
		this.x = ran.nextInt((int) maxX);
		this.y = ran.nextInt((int) maxY);
	}
	public Coordinate(float x, float y) {
		this.x = x;
		this.y = y;
	}
	public float distanceFrom(Coordinate c) {
		return (float)(Math.sqrt( Math.pow((this.y - c.y), 2) + Math.pow(this.x - c.x, 2)));
	}
	public void print() {
		System.out.println("(" + this.x + "," + this.y + ")");
	}
}
class Vertices {
	private Vertex[] vertices;
	public Vertices(int vertexCount) {
		vertices = new Vertex[vertexCount];
		for(int i = 0; i < vertexCount; i++) {
			vertices[i] = (new Vertex(i));
		}
	}
	
	public int bruteForceAlgorithmCost(int[][] costMatrix) {
		int cheapestCost = -1;
		ArrayList<ArrayList<Integer>> permutations = this.permutate();
		for(int i = 0; i < permutations.size(); i++) {
			//If the first item is not zero, then skip..
			if(permutations.get(i).get(0) != 0) {
				continue;
			} else {
				//Calculate the cost of the permutation
				int cost = 0; 
				for(int j = 0; j < permutations.get(i).size(); j++) {
					if(permutations.get(i).get(j)+1 >= permutations.get(i).size())
						cost += costMatrix[permutations.get(i).get(j)][0];
					else
						cost += costMatrix[permutations.get(i).get(j)][permutations.get(i).get(j)+1]; 
				}
				if(cheapestCost == -1 || cheapestCost > cost) {
					cheapestCost = cost;
				}
			}
		}
		return cheapestCost;		
	}
	
	public int greedyAlgorithmCost(int[][] costMatrix) {
		//From 0 and back..
		int lastVertex = 0;
		int vertex = 0;
		int cost = 0;
		while(vertex != -1) {
			lastVertex = vertex;
			vertex = findCheapestUnvisitedVertex(costMatrix, vertex);
			if(vertex != -1) {
				cost+=costMatrix[lastVertex][vertex];
				//System.out.println(lastVertex + "=>" + vertex + ":" + costMatrix[lastVertex][vertex]);
				this.vertices[vertex].visit();
			}
		}
		//Find the cost from last vertex back to 0...
		cost += costMatrix[lastVertex][0];
		//System.out.println(lastVertex + "=>" + 0 + ":" + costMatrix[lastVertex][0]);

		return cost;
	}
	
	/**
	 * 
	 * @param costMatrix
	 * @param a
	 * @return The index of the vertex for the cheapest cost from point a.. -1 if all are visited.
	 */
	private int findCheapestUnvisitedVertex(int[][] costMatrix, int a) {
		int cheapestIndex = -1;
		for(int i = 1; i < this.vertices.length; i++) {
			if((cheapestIndex == -1) || (costMatrix[a][i] < costMatrix[a][cheapestIndex])) {
				if(!this.vertices[i].visited)
					cheapestIndex = this.vertices[i].index;
			}
		}
		return cheapestIndex;
	}
	
	public ArrayList<ArrayList<Integer>> permutate() {		
		int[] input = new int[vertices.length];
		for(int i = 0; i < this.vertices.length; i++) {
			input[i] = this.vertices[i].index;
		}
        int n = input.length; 
        Permutation permutation = new Permutation(); 
        return permutation.permute(input, 0, n-1); 
/*        
        for(int i = 0; i < permutation.permutations.size(); i++) {
        	for(int j = 0; j < permutation.permutations.get(i).size(); j++) {
        		//System.out.print(permutation.permutations.get(i).get(j) + ",");
        	}
        	//System.out.println();
        }*/
	}
}
class Vertex {
	public final int index;
	public boolean visited = false;
	public Vertex(int index) {
		this.index = index;
	}
	public boolean visit() {
		return this.visited = true;
	}
}
class Permutation 
{
	public ArrayList<ArrayList<Integer>> permutations = new ArrayList<ArrayList<Integer>>();
    /** 
    * permutation function 
    * @param str string to calculate permutation for 
    * @param l starting index 
    * @param r end index 
    */
    ArrayList<ArrayList<Integer>> permute(int[] str, int l, int r) 
    { 
        if (l == r) {
        	ArrayList<Integer> list = new ArrayList<Integer>();
            for(int i = 0; i < str.length; i++) {
            	list.add(str[i]);
            }
            this.permutations.add(list);
        }
        else
        { 
            for (int i = l; i <= r; i++) 
            { 
                str = swap(str,l,i); 
                permute(str, l+1, r); 
                str = swap(str,l,i); 
            } 
        } 
        return this.permutations;
    } 
 
    /** 
    * Swap Characters at position 
    * @param a string value 
    * @param i position 1 
    * @param j position 2 
    * @return swapped string 
    */
    public int[] swap(int[] a, int i, int j) 
    { 
        int temp; 
        int[] charArray = a; 
        temp = charArray[i] ; 
        charArray[i] = charArray[j]; 
        charArray[j] = temp; 
        return charArray; 
    } 
 
} 