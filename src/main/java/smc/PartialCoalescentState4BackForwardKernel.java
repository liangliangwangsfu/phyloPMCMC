package smc;

import pty.io.Dataset;
import pty.smc.PartialCoalescentState;
import pty.smc.models.CTMC;

public class PartialCoalescentState4BackForwardKernel extends PartialCoalescentState {
	private PartialCoalescentState4BackForwardKernel parent = null;
	private double deltaOld = 0;
	private int[] oldIndx=new int[2];
	private int[] newIndx1=new int[2];
	private int[] newIndx2=new int[2];
	
	public double getDeltaOld() {
		return deltaOld;
	}

	public void setDeltaOld(double deltaOld) {
		this.deltaOld = deltaOld;
	}

	public PartialCoalescentState4BackForwardKernel(PartialCoalescentState pcs,
			PartialCoalescentState4BackForwardKernel parent, double deltaOld, int[] oldIndx, int[] newIndx1, int[] newIndx2) {
		super(pcs);
		this.parent = parent;
		this.deltaOld = deltaOld;
		this.setOldIndx(oldIndx);
		this.setNewIndx1(newIndx1);
		this.setNewIndx1(newIndx2);
	}
	
	public static PartialCoalescentState initFastState(Dataset data, CTMC ctmc,
			boolean isClock) {
		PartialCoalescentState.initFastState(false, data, ctmc, isClock);
		return initFastState(false, data, ctmc, isClock);
	}

	public PartialCoalescentState4BackForwardKernel parentState() {
		return parent;
	}

	public void setParent(PartialCoalescentState4BackForwardKernel parent) {
		this.parent = parent;
	}

	public boolean hasParent()
{
		return (this.parent != null);
			}

	public int[] getOldIndx() {
		return oldIndx;
	}

	public int getOldIndxLeft() {
		return oldIndx[0];
	}

	public int getOldIndxRight() {
		return oldIndx[1];
	}
	
	
	public void setOldIndx(int[] oldIndx) {
		this.oldIndx = oldIndx;
	}

	
	public int[] getNewIndx1() {
		return newIndx1;
	}

	public int getNewIndx1Left() {
		return newIndx1[0];
	}

	public int getNewIndx1Right() {
		return newIndx1[1];
	}
	
	public void setNewIndx1(int[] newIndx1) {
		this.newIndx1 = newIndx1;
	}

	public int[] getNewIndx2() {
		return newIndx2;
	}
	
	public int getNewIndx2Left() {
		return newIndx1[0];
	}

	public int getNewIndx2Right() {
		return newIndx1[1];
	}

	public void setNewIndx2(int[] newIndx2) {
		this.newIndx2 = newIndx2;
	}

}
