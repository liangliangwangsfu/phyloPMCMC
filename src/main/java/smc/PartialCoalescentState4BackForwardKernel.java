package smc;

import static nuts.util.CollUtils.list;
import java.util.List;
import fig.basic.Pair;
import pty.io.Dataset;
import pty.smc.PartialCoalescentState;
import pty.smc.models.CTMC;

public class PartialCoalescentState4BackForwardKernel{
	private  PartialCoalescentState currentState=null;
	private  PartialCoalescentState previousState=null;
	private  PartialCoalescentState midState=null;
	
	public PartialCoalescentState getCurrentState() {
		return currentState;
	}

	public void setCurrentState(PartialCoalescentState currentState) {
		this.currentState = currentState;
	}

	public PartialCoalescentState getPreviousState() {
		return previousState;
	}
	
	public void setPreviousState(PartialCoalescentState previousState) {
		this.previousState = previousState;
	}
	
	public PartialCoalescentState getMidState() {
		return midState;
	}

	public void setMidState(PartialCoalescentState midState) {
		this.midState = midState;
	}

	private PartialCoalescentState4BackForwardKernel parent = null;
	private double delta = 0;
	private int[] indxState=new int[2];
		

	public PartialCoalescentState4BackForwardKernel(PartialCoalescentState currentState,
			PartialCoalescentState midState, PartialCoalescentState previousState,			
			PartialCoalescentState4BackForwardKernel parent, double delta, int[] indxState) {
		this.currentState=currentState;
		this.midState=midState;
		this.previousState=previousState;
		this.parent = parent;
		this.setDelta(delta);
		this.setIndxState(indxState);
	}

	public static PartialCoalescentState initFastState(Dataset data, CTMC ctmc,
			boolean isClock) {
		return PartialCoalescentState.initFastState(false, data, ctmc, isClock);		
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
	
	
	public static double forwardDensity(PartialCoalescentState4BackForwardKernel thisState, PartialCoalescentState4BackForwardKernel newState)
	{
		double result=Double.NEGATIVE_INFINITY;		
		if(thisState.getMidState() == newState.getPreviousState())
		{
			result=newState.getCurrentState().logLikelihoodRatio()
					+ newState.getMidState().logLikelihoodRatio()
					- thisState.getCurrentState().logLikelihoodRatio();	        
		}
		return result;		
	}
	
	public static List<Pair<PartialCoalescentState4BackForwardKernel,Double>>  restoreSequence(PartialCoalescentState4BackForwardKernel finalState)
	{	
		List<Pair<PartialCoalescentState4BackForwardKernel,Double>> result=list();		
		PartialCoalescentState4BackForwardKernel current=finalState;				
		while(current.parentState().hasParent())
		{
			result.add(Pair.makePair(current, current.getCurrentState().logLikelihoodRatio() + current.getMidState().logLikelihoodRatio()-  current.parentState().getCurrentState().logLikelihoodRatio())); 
		    current=current.parentState();		
		}
		result.add(Pair.makePair(current, current.getCurrentState().logLikelihoodRatio()));
		List<Pair<PartialCoalescentState4BackForwardKernel,Double>> finalResult=list();
		for(int i=result.size()-1;i>=0;i--)
			finalResult.add(result.get(i));
		return finalResult;
	}
	
	
	public static void main(String [] args)
	  {
//        File data = new File("/Users/liangliangwang/eclipse-workspace/phyloPMCMC/state/execs/7.exec/output/sim-379066948.msf");	
//		Dataset dataset = DatasetUtils.fromAlignment(data, SequenceType.DNA);
//		CTMC ctmc = CTMC.SimpleCTMC.dnaCTMC(dataset.nSites());
//		PartialCoalescentState init0 = PartialCoalescentState
//				.initFastState(dataset, ctmc, true);
//		PartialCoalescentState4BackForwardKernel init = new PartialCoalescentState4BackForwardKernel(
//				init0, null, 0, new int[] {-1,-1});
//
//		ParticleKernel<PartialCoalescentState4BackForwardKernel> kernel = new BackForwardKernel(
//				init);
//
//		ParticleFilter<PartialCoalescentState4BackForwardKernel> pf = new ParticleFilter<PartialCoalescentState4BackForwardKernel>();
//		pf.nThreads = 1;
//		pf.resampleLastRound = false;
//		pf.N=100;
//		StoreProcessor<PartialCoalescentState4BackForwardKernel> pro = new StoreProcessor<PartialCoalescentState4BackForwardKernel>();
//		pf.sample(kernel, pro);
//		PartialCoalescentState4BackForwardKernel sampled = pro.sample(new Random(441));		
//		List<Pair<PartialCoalescentState4BackForwardKernel,Double>> conditionedPath = restoreSequence(sampled);
	  }


	public double getDelta() {
		return delta;
	}


	public void setDelta(double delta) {
		this.delta = delta;
	}


	public int[] getIndxState() {
		return indxState;
	}


	public void setIndxState(int[] indxState) {
		this.indxState = indxState;
	}

}
