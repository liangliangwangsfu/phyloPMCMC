package smc;
import java.util.List;
import java.util.Random;

import nuts.math.Sampling;
import pty.smc.LazyParticleFilter.LazyParticleKernel;
import pty.smc.PartialCoalescentState;
import pty.smc.ParticleKernel;
import fig.basic.Option;
import fig.basic.Pair;

/**
 * 
 * @author Liangliang Wang
 */
public class BackForwardKernel2 implements
LazyParticleKernel<PartialCoalescentState4BackForwardKernel2>,ParticleKernel<PartialCoalescentState4BackForwardKernel2> {
	@Option
	public static boolean printBranchLengthMagnitudes = false;
	private final PartialCoalescentState4BackForwardKernel2 initial;

	public BackForwardKernel2(PartialCoalescentState4BackForwardKernel2 initial) {
		this.initial = initial;
	}

	public PartialCoalescentState4BackForwardKernel2 getInitial() {
		return initial;
	}
		

	public int nIterationsLeft(
			PartialCoalescentState4BackForwardKernel2 partialState) {
	//	System.out.println(partialState.getCurrentState().toString());
		return partialState.getCurrentState().nIterationsLeft();
	}

	public Object _next0(Random rand,
			PartialCoalescentState4BackForwardKernel2 current,
			boolean isPeek) {

		final double delta = Sampling.sampleExponential(rand,
				(initial.getCurrentState().isClock() ? 0.1 : 0.05) / nChoose2(current.getCurrentState().nRoots()));

		// 2- sample a random pair (without replacement)
		List<Integer> sampledIndices = Sampling.sampleWithoutReplacement(rand,
				current.getCurrentState().nRoots(), 2);

		final int i0 = sampledIndices.get(0), i1 = sampledIndices.get(1);

		double leftIncrement = 0.0, rightIncrement = 0.0;
		if (!initial.getCurrentState().isClock()) {
			double incr = Sampling.sampleExponential(rand, 0.1);
			if (rand.nextBoolean())
				leftIncrement += incr;
			else
				rightIncrement += incr;
		}

		PartialCoalescentState4BackForwardKernel2 result = null;
		Double logw = null;
		if (isPeek)
			logw = current.getCurrentState().peekLogLikelihoodRatio(i0, i1, delta, leftIncrement,
					rightIncrement);
		else
		{
			PartialCoalescentState resultpcs = current.getCurrentState()
					.coalesce(i0, i1, delta, leftIncrement,
							rightIncrement);
			result = new PartialCoalescentState4BackForwardKernel2(resultpcs,current.getCurrentState(),null,
					current, delta, new int[]{i0, i1});
		}

		if (isPeek)
			return logw;
		else
			return Pair.makePair(result, result.getCurrentState().logLikelihoodRatio());
	}

	public Object _next1(Random rand,
			PartialCoalescentState4BackForwardKernel2 current,
			boolean isPeek) {
		boolean self=false;		
		double oldLogLikelihoodRatio = 0;
		double logExpDensityDeltaOld = 0;
		// 1- go back to current state's parent partial state
		PartialCoalescentState4BackForwardKernel2 parent = current.parentState();		
		//PartialCoalescentState4BackForwardKernel2 result0 = null;
		PartialCoalescentState result0 = null;
		
		// 3- sample a random pair (without replacement)
		List<Integer> sampledIndices0 = Sampling.sampleWithoutReplacement(rand,
				parent.getCurrentState().nRoots(), 2);
		final int i00 = sampledIndices0.get(0), i01 = sampledIndices0.get(1);
		int[] currentIndx=current.getIndxState();
		if((i00==currentIndx[0] && i01==currentIndx[1]) || (i00==currentIndx[1] && i01==currentIndx[0] ))
		{
			self=true;
			result0=current.getCurrentState();
		}
		else {
//		System.out.println(parent.toString());;
		oldLogLikelihoodRatio = current.getCurrentState().logLikelihoodRatio();
		double deltaOld = current.getDelta();

		// 2- sample the exponential waiting time
		
		double param0=(initial.getCurrentState().isClock() ? 0.1 : 0.05)
				/ nChoose2(parent.getCurrentState().nRoots());
		final double delta0 = Sampling.sampleExponential(rand, param0);
		double leftIncrement0 = 0.0, rightIncrement0 = 0.0;
		//		if (!initial.isClock()) {
		//			double incr0 = Sampling.sampleExponential(rand, 0.1), incr1 = Sampling
		//					.sampleExponential(rand, 0.1);
		//			if (rand.nextBoolean()) {
		//				leftIncrement0 += incr0;
		//			}
		// else {
		//				rightIncrement0 += incr0;
		//			}
		//		}
		logExpDensityDeltaOld = Sampling.exponentialLogDensity(param0,deltaOld);				
		
		//PartialCoalescentState resultMid
		result0=current.getMidState().coalesce(i00, i01, delta0, leftIncrement0, rightIncrement0); 
		
		//result = new PartialCoalescentState4BackForwardKernel2(resultMid.coalesce(i10, i11, delta1, leftIncrement1, rightIncrement1),result0.getCurrentState(),result0.getMidState(),result0,delta1,new int[]{i10, i11});
				
	//	result0 = new PartialCoalescentState4BackForwardKernel2(current.getPreviousState().coalesce(i00, i01, delta0, leftIncrement0, rightIncrement0), parent.getCurrentState(),parent.getMidState(), parent, delta0, new int[]{i00, i01});
	//	result0 = new PartialCoalescentState4BackForwardKernel2(current.getPreviousState().coalesce(i00, i01, delta0, leftIncrement0, rightIncrement0),resultMid,parent.getMidState(), parent, delta0, new int[]{i00, i01});
		}
		List<Integer> sampledIndices1 = Sampling.sampleWithoutReplacement(rand,
				result0.nRoots(), 2);
		final int i10 = sampledIndices1.get(0), i11 = sampledIndices1.get(1);
		double param1 = (initial.getCurrentState().isClock() ? 0.1 : 0.05)
				/ nChoose2(current.getCurrentState().nRoots());
		final double delta1 = Sampling.sampleExponential(rand, param1);		
		double leftIncrement1 = 0.0, rightIncrement1 = 0.0;
		PartialCoalescentState4BackForwardKernel2 result = null;
		Double loglikeRatio0=result0.logLikelihoodRatio();
		Double logw =null;
		if (isPeek) {
			if(self)
				logw= result0.peekLogLikelihoodRatio(i10, i11, delta1,
						leftIncrement1, rightIncrement1);
			else
				logw= loglikeRatio0+result0.peekLogLikelihoodRatio(i10, i11, delta1,
					leftIncrement1, rightIncrement1) - oldLogLikelihoodRatio - logExpDensityDeltaOld ;
		} else {		
			result = new PartialCoalescentState4BackForwardKernel2(result0.coalesce(i10, i11, delta1, leftIncrement1, rightIncrement1),result0,current.getMidState(),current,delta1,new int[]{i10, i11});
		}

		// 4- the weight update is simply equal to the ratio of the new
		// likelihood score to the old one
		if (isPeek)
			return logw;
		else{
			if(self)return Pair.makePair(result, result.getCurrentState().logLikelihoodRatio());
			else
			return Pair.makePair(result, result.getCurrentState().logLikelihoodRatio()
					+ loglikeRatio0- oldLogLikelihoodRatio - logExpDensityDeltaOld);
		}
	}

	public static double nChoose2(double n) {
		return n * (n - 1) / 2;
	}


	public Object _next(Random rand,
			PartialCoalescentState4BackForwardKernel2 current, boolean isPeek) {
		//System.out.println(current.hasParent());
		if (current.hasParent())
			return _next1(rand, current, isPeek);
		else
			return _next0(rand, current, isPeek);
	}

	@Override
	public Pair<PartialCoalescentState4BackForwardKernel2, Double> next(
			Random rand, PartialCoalescentState4BackForwardKernel2 current) {
		return (Pair) _next(rand, current, false);

	}

	@Override
	public double peekNext(Random rand,
			PartialCoalescentState4BackForwardKernel2 current) {
		return (Double) _next(rand, current, true);
	}	

}
