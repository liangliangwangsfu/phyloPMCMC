package smc;

import static nuts.util.CollUtils.list;

import java.io.File;
import java.util.List;
import java.util.Random;

import fig.basic.Pair;
import ma.SequenceType;
import nuts.math.Sampling;
import pty.io.Dataset;
import pty.io.Dataset.DatasetUtils;
import pty.smc.PartialCoalescentState;
import pty.smc.ParticleFilter;
import pty.smc.ParticleFilter.StoreProcessor;
import pty.smc.ParticleKernel;
import pty.smc.models.CTMC;

public class PartialCoalescentState4BackForwardKernel extends PartialCoalescentState {
	private PartialCoalescentState4BackForwardKernel parent = null;
	private double latestDelta = 0;
	private int[] indxInParentState=new int[2];

	public PartialCoalescentState4BackForwardKernel(PartialCoalescentState pcs,
			PartialCoalescentState4BackForwardKernel parent, double latestDelta, int[] indxInParentState) {
		super(pcs);
		this.parent = parent;
		this.setLatestDelta(latestDelta);
		this.setIndxInParentState(indxInParentState);
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



	public double getLatestDelta() {
		return latestDelta;
	}

	public void setLatestDelta(double latestDelta) {
		this.latestDelta = latestDelta;
	}

	public int[] getIndxInParentState() {
		return indxInParentState;
	}

	public void setIndxInParentState(int[] indxInParentState) {
		this.indxInParentState = indxInParentState;
	}
	
	
	public static double forwardDensity(PartialCoalescentState4BackForwardKernel thisState, PartialCoalescentState4BackForwardKernel newState)
	{
		double result=0;		
		PartialCoalescentState4BackForwardKernel grandmaOfNewState=newState.parentState().parentState();
		PartialCoalescentState4BackForwardKernel parentOfthisState=thisState.parentState();
		if(grandmaOfNewState.equals(parentOfthisState))
		{
			double numRootPairs0=BackForwardKernel.nChoose2(grandmaOfNewState.nRoots());
			double param0= 0.1 / numRootPairs0;				
	        double logExpDensityDelta0 = Sampling.exponentialLogDensity(param0, newState.parentState().getLatestDelta());	        
			double numRootPairs1=BackForwardKernel.nChoose2(newState.parentState().nRoots());
			double param1= 0.1 / numRootPairs1;				
	        double logExpDensityDelta1 = Sampling.exponentialLogDensity(param1, newState.getLatestDelta());	        
	        result=logExpDensityDelta0+logExpDensityDelta1-Math.log(numRootPairs0)-Math.log(numRootPairs1);
		}		
		return result;		
	}
	
	public static List<Pair<PartialCoalescentState4BackForwardKernel,Double>>  restoreSequence(	PartialCoalescentState4BackForwardKernel finalState)
	{
	
		List<Pair<PartialCoalescentState4BackForwardKernel,Double>> result=list();		
		PartialCoalescentState4BackForwardKernel current=finalState;				
		while(current.parentState().hasParent())
		{
			double param0 = 0.1 / BackForwardKernel.nChoose2(current.parentState().parentState().nRoots());
			final double delta0 = current.parentState().parentState().getLatestDelta();			
			double logExpDensityDeltaOld = Sampling.exponentialLogDensity(param0,
					delta0);		
			result.add(Pair.makePair(current, current.logLikelihoodRatio() + current.parentState().logLikelihoodRatio()-  current.parentState().parentState().logLikelihoodRatio() - logExpDensityDeltaOld));
//		System.out.println(i++);
//		System.out.println(current.parentState().toString());
//		System.out.println(current.getLatestDelta()+": "+current.getIndxInParentState()[0]+","+current.getIndxInParentState()[1]);
//		System.out.println(current.toString());
		current=current.parentState();		
		}
		result.add(Pair.makePair(current, current.logLikelihoodRatio()));
		//System.out.println(result.size());
		
		List<Pair<PartialCoalescentState4BackForwardKernel,Double>> finalResult=list();
		for(int i=result.size()-1;i>=0;i--)
			finalResult.add(result.get(i));
		return finalResult;
	}
	
	
	
	
	
	public static void main(String [] args)
	  {
        File data = new File("/Users/liangliangwang/eclipse-workspace/phyloPMCMC/state/execs/7.exec/output/sim-379066948.msf");	
		Dataset dataset = DatasetUtils.fromAlignment(data, SequenceType.DNA);
		CTMC ctmc = CTMC.SimpleCTMC.dnaCTMC(dataset.nSites());
		PartialCoalescentState init0 = PartialCoalescentState
				.initFastState(dataset, ctmc, true);
		PartialCoalescentState4BackForwardKernel init = new PartialCoalescentState4BackForwardKernel(
				init0, null, 0, new int[] {-1,-1});

		ParticleKernel<PartialCoalescentState4BackForwardKernel> kernel = new BackForwardKernel(
				init);

		ParticleFilter<PartialCoalescentState4BackForwardKernel> pf = new ParticleFilter<PartialCoalescentState4BackForwardKernel>();
		pf.nThreads = 1;
		pf.resampleLastRound = false;
		pf.N=100;
		StoreProcessor<PartialCoalescentState4BackForwardKernel> pro = new StoreProcessor<PartialCoalescentState4BackForwardKernel>();
		pf.sample(kernel, pro);
		PartialCoalescentState4BackForwardKernel sampled = pro.sample(new Random(441));		
		List<Pair<PartialCoalescentState4BackForwardKernel,Double>> conditionedPath = restoreSequence(sampled);
	  }

}
