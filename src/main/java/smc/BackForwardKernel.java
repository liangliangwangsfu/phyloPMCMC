package smc;
import java.util.List;
import java.util.Random;

import nuts.math.Sampling;
import nuts.util.Arbre;
import pty.smc.LazyParticleFilter.LazyParticleKernel;
import pty.smc.PartialCoalescentState.CoalescentNode;
import pty.smc.PartialCoalescentState;
import pty.smc.ParticleKernel;
import fig.basic.Option;
import fig.basic.Pair;

/**
 * 
 * @author Liangliang Wang
 */
public class BackForwardKernel implements
LazyParticleKernel<PartialCoalescentState4BackForwardKernel>,ParticleKernel<PartialCoalescentState4BackForwardKernel> {
	@Option
	public static boolean printBranchLengthMagnitudes = false;
	private final PartialCoalescentState4BackForwardKernel initial;
	

	public BackForwardKernel(PartialCoalescentState4BackForwardKernel initial) {
		this.initial = initial;
	}

	public PartialCoalescentState4BackForwardKernel getInitial() {
		return initial;
	}


	public int nIterationsLeft(
			PartialCoalescentState4BackForwardKernel partialState) {
		//	System.out.println(partialState.getCurrentState().toString());
		return partialState.getCurrentState().nIterationsLeft();
	}
	


	public Object _next0(Random rand,
			PartialCoalescentState4BackForwardKernel current,
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

		PartialCoalescentState4BackForwardKernel result = null;
		Double logw = null;
		if (isPeek)
			logw = current.getCurrentState().peekLogLikelihoodRatio(i0, i1, delta, leftIncrement,
					rightIncrement);
		else
		{
			PartialCoalescentState resultpcs = current.getCurrentState()
					.coalesce(i0, i1, delta, leftIncrement,
							rightIncrement);
			result = new PartialCoalescentState4BackForwardKernel(resultpcs,current.getCurrentState(),null,
					current, delta, new int[]{i0, i1});
		}

		if (isPeek)
			return logw;
		else
			return Pair.makePair(result, result.getCurrentState().logLikelihoodRatio());
	}

	public Object _next1(Random rand,
			PartialCoalescentState4BackForwardKernel current,
			boolean isPeek) {
		//boolean self=false;		
		double oldLogLikelihoodRatio = 0;
		double logExpDensityDeltaOld = 0;
		// 1- go back to current state's parent partial state
		PartialCoalescentState4BackForwardKernel parent = current.parentState();		
		//PartialCoalescentState4BackForwardKernel2 result0 = null;
		PartialCoalescentState result00 = current.getCurrentState();
		PartialCoalescentState result0 = null;
		PartialCoalescentState result1 = null;

		// 3- sample a random pair (without replacement)
		List<Integer> sampledIndices0 = Sampling.sampleWithoutReplacement(rand,
				parent.getCurrentState().nRoots(), 2);
		final int i00 = sampledIndices0.get(0), i01 = sampledIndices0.get(1);
		int[] currentIndx=current.getIndxState();
		
			//		System.out.println(parent.toString());;
			oldLogLikelihoodRatio = current.getCurrentState().logLikelihoodRatio();
			double deltaOld = current.getDelta();
			

			// 2- sample the exponential waiting time

			double param0=(initial.getCurrentState().isClock() ? 0.1 : 0.05)
					/ nChoose2(parent.getCurrentState().nRoots());
			double tempdelta0 = Sampling.sampleExponential(rand, param0);
			int deltaIter = 0;
			while(tempdelta0 < Math.pow(10, -10) && deltaIter < 5) {
				tempdelta0 = Sampling.sampleExponential(rand, param0);
				deltaIter = deltaIter +1;				
			}
			
			if(tempdelta0 < Math.pow(10, -10)) {
				tempdelta0 = Math.pow(10, -10);
			}
			
			final double delta0 = tempdelta0;
			
			double leftIncrement0 = 0.0, rightIncrement0 = 0.0;


			result0=current.getMidState().coalesce(i00, i01, delta0, leftIncrement0, rightIncrement0); 

		List<Integer> sampledIndices1 = Sampling.sampleWithoutReplacement(rand,
				result0.nRoots(), 2);
		final int i10 = sampledIndices1.get(0), i11 = sampledIndices1.get(1);
		final double nPossiblePairs = nChoose2(current.getCurrentState().nRoots());
		double param1 = (initial.getCurrentState().isClock() ? 0.1 : 0.05)/nPossiblePairs;
		final double delta1 = Sampling.sampleExponential(rand, param1);
		double leftIncrement1 = 0.0, rightIncrement1 = 0.0;
		PartialCoalescentState4BackForwardKernel result = null;
		Double loglikeRatio0=result0.logLikelihoodRatio();
		result1 = result0.coalesce(i10, i11, delta1, leftIncrement1, rightIncrement1); 
		
		Double logw =null;

		logw = loglikeRatio0 + result1.logLikelihoodRatio() - oldLogLikelihoodRatio;

		result = new PartialCoalescentState4BackForwardKernel(result1,result0,current.getMidState(),current,delta1,new int[]{i10, i11});       

		if (isPeek)
			return logw; 
		else{
			return Pair.makePair(result, result.getCurrentState().logLikelihoodRatio() 
								+ loglikeRatio0 - oldLogLikelihoodRatio);
			}
		}

//	}

	public static double nChoose2(double n) {
		return n * (n - 1) / 2;
	}


	public Object _next(Random rand,
			PartialCoalescentState4BackForwardKernel current, boolean isPeek) {
		//System.out.println(current.hasParent());
		if (current.hasParent())
			return _next1(rand, current, isPeek);
		else
			return _next0(rand, current, isPeek);
	}

	@Override
	public Pair<PartialCoalescentState4BackForwardKernel, Double> next(
			Random rand, PartialCoalescentState4BackForwardKernel current) {
		return (Pair) _next(rand, current, false);

	}

	@Override
	public double peekNext(Random rand,
			PartialCoalescentState4BackForwardKernel current) {
		return (Double) _next(rand, current, true);
	}	

}
