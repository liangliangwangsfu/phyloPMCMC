package smc;

import java.util.*;
import pty.smc.PartialCoalescentState;
import pty.smc.ParticleKernel;
import monaco.process.ProcessSchedule;
import monaco.process.ProcessScheduleContext;
import monaco.process.ResampleStatus;
import nuts.lang.ArrayUtils;
import nuts.math.Fct;
import nuts.math.Id;
import nuts.math.MeasureZeroException;
import nuts.math.Sampling;
import nuts.maxent.SloppyMath;
import nuts.util.CollUtils;
import nuts.util.Counter;
import nuts.util.Hasher;
import fig.basic.LogInfo;
import fig.basic.NumUtils;
import fig.basic.Option;
import fig.basic.Pair;
import fig.basic.Parallelizer;
import fig.prob.Multinomial;
import fig.prob.SampleUtils;
import goblin.BayesRiskMinimizer;
import goblin.BayesRiskMinimizer.LossFct;

public final class PGASParticleFilter<S> {
	@Option
	public boolean verbose = false;
	@Option
	public int N = 100;
	@Option
	public Random rand = new Random(1);
	@Option
	public boolean resampleLastRound = true;
	@Option
	public int nThreads = 1;
	@Option
	public ResamplingStrategy resamplingStrategy = ResamplingStrategy.ALWAYS;
	public static final double essRatioThreshold = 0.5;

	@Option
	public boolean usePGAS = true;

	private List<S> conditional = null; // exclude initial state
	private double[] conditionalUnnormWeights = null;

	private ProcessSchedule schedule = null;

	public void setProcessSchedule(ProcessSchedule schedule) {
		this.schedule = schedule;
	}

	public void setConditional(List<S> conditional, double[] conditionalUnnormWeights) {
		// no I think it's ok
		// if (true)
		// throw new RuntimeException("seems like conditional not used in
		// resampling step?");
		if (conditional.size() != conditionalUnnormWeights.length)
			throw new RuntimeException();
		this.conditional = conditional;
		this.conditionalUnnormWeights = conditionalUnnormWeights;
	}

	public void setUnconditional() {
		this.conditional = null;
		this.conditionalUnnormWeights = null;
	}

	public boolean isConditional() {
		return conditional != null;
	}

	public static enum ResamplingStrategy {
		ALWAYS {
			@Override
			public boolean needResample(double[] w) {
				return true;
			}
		},
		NEVER {
			@Override
			public boolean needResample(double[] w) {
				return false;
			}
		},
		ESS {
			@Override
			public boolean needResample(final double[] weights) {
				final double ess = ess(weights);
				final double threshold = essRatioThreshold * weights.length;
				return ess < threshold;
			}
		};
		abstract boolean needResample(double[] weights);
	}

	public static double ess(double[] ws) {
		double sumOfSqr = 0.0;
		for (double w : ws)
			sumOfSqr += w * w;
		return 1.0 / sumOfSqr;
	}

	public static <S> void bootstrapFilter(final ParticleKernel<S> kernel, final ParticleProcessor<S> processor,
			final int N, final Random rand) {
		PGASParticleFilter<S> PF = new PGASParticleFilter<S>();
		PF.N = N;
		PF.rand = rand;
		PF.resampleLastRound = false;
		PF.resamplingStrategy = ResamplingStrategy.ALWAYS;
		PF.sample(kernel, processor);
	}

	private long[] seeds;
	private List<S> samples;
	private double[] logWeights;
	private double[] normalizedWeights;
	// private List<Integer> ancestors;

	public List<S> getSamples() {
		return samples;
	}

	public double[] getLogWeights() {
		return logWeights;
	}


	private void propagateAndComputeWeights(final ParticleKernel<S> kernel, final int t) {
		if (verbose)
			LogInfo.track("Processing...", false);
		// System.out.println("Creating seed=" + rand.nextDouble());
		seeds = Sampling.createSeeds(N, rand); // so that result for a given
		// randomization is
		// #-of-threads-invariant
		Parallelizer<Integer> parallelizer = new Parallelizer<Integer>(nThreads);
		parallelizer.setPrimaryThread();
		parallelizer.process(CollUtils.ints(N), new Parallelizer.Processor<Integer>() {
			public void process(Integer x, int _i, int _n, boolean log) {
				if (log && ((x + 1) % 5 == 0))
					if (verbose)
						LogInfo.logs("Particle " + (x + 1) + "/" + N);
				Random rand = new Random(seeds[x]);
				if (x == 0 && isConditional()) {
					samples.set(x, conditional.get(t));
					logWeights[x] = conditionalUnnormWeights[t];
				} else {
					final Pair<S, Double> current = kernel.next(rand, samples.get(x));
					if (current == null) {
						samples.set(x, null);
						logWeights[x] = Double.NEGATIVE_INFINITY;
					} else {
						samples.set(x, current.getFirst());
						logWeights[x] = Math.log(normalizedWeights[x]) + current.getSecond();
					}
				}
			}
		});
		if (verbose)
			LogInfo.end_track();
	}

	private double lognorm = 0.0;

	public double estimateNormalizer() {
		// if (resamplingStrategy != ResamplingStrategy.ALWAYS)
		// throw new RuntimeException();
		return lognorm;
	}

	private void init(final ParticleKernel<S> kernel) {
		// ancestors = new ArrayList<Integer>(N);
		// for (int n = 0; n < N; n++)
		// ancestors.add(n);
		lognorm = 0.0;
		samples = new ArrayList<S>(N);
		S initial = kernel.getInitial();
		for (int n = 0; n < N; n++)
			samples.add(initial);
		logWeights = new double[N]; // init to zero
		normalizedWeights = new double[N];
		for(int i=0;i<normalizedWeights.length;i++) normalizedWeights[i]=1.0/N;
	}

	private void newProcess(int t, double[] normalizedWeights, ParticleProcessor<S> processor, int T) {
		if ((schedule == null && t == T - 1) || (schedule != null
				&& schedule.shouldProcess(new ProcessScheduleContext(t, t == T - 1, ResampleStatus.NA)))) {
			if (verbose)
				LogInfo.track("Processing particles");
			for (int n = 0; n < normalizedWeights.length; n++) {
				if (verbose)
					LogInfo.logs("Particle " + (n + 1) + "/" + normalizedWeights.length);
				if (samples.get(n) != null) {
					processor.process(samples.get(n), normalizedWeights[n]);
					// samples.set(n, null);
				}
			}
			if (verbose)
				LogInfo.end_track();
		}
	}

	/**
	 * @param <S>
	 * @param kernel
	 *            model
	 * @param processor
	 *            what to do with the produced sample
	 */
	public void sample(final ParticleKernel<S> kernel, final ParticleProcessor<S> processor) {
		init(kernel);		

		final int T = kernel.nIterationsLeft(kernel.getInitial());

		// System.out.println("T=" + T);
		if (isConditional() && conditional.size() != T)
			throw new RuntimeException();
		for (int t = 0; t < T; t++) {
			if (verbose)
				LogInfo.track("Particle generation " + (t + 1) + "/" + T, true);
			propagateAndComputeWeights(kernel, t);
			if(t>0)lognorm += SloppyMath.logAdd(logWeights);
			normalizedWeights = logWeights.clone();
			NumUtils.expNormalize(normalizedWeights);
			if (verbose)
				LogInfo.logs("LargestNormalizedWeights=" + ArrayUtils.max(normalizedWeights));
			if (verbose)
				LogInfo.logs("RelativeESS=" + ess(normalizedWeights) / normalizedWeights.length);
			// if (verbose) LogInfo.logs("NumberOfAncestors=" + new
			// HashSet<Integer>(ancestors).size());			
			newProcess(t, normalizedWeights, processor, T);
			if (t < T - 1 && (hasNulls(samples) || resamplingStrategy.needResample(normalizedWeights))) {

				if (usePGAS && isConditional()
						&& conditional.get(t) instanceof PartialCoalescentState4BackForwardKernel2) {					
					double[] forwardDensityWeights = new double[N];
					PartialCoalescentState4BackForwardKernel2 conditionedState = (PartialCoalescentState4BackForwardKernel2) conditional
							.get(t+1);
					for (int k = 0; k < N; k++) {						 
						forwardDensityWeights[k]= PartialCoalescentState4BackForwardKernel2.forwardDensity(
								(PartialCoalescentState4BackForwardKernel2) samples.get(k), conditionedState);
						//	if(tmp!=0) forwardDensityWeights[k]=tmp+logWeights[k];													
						//		if(forwardDensityWeights[k]!=0)System.out.print(k+": "+forwardDensityWeights[k]+"	"+tmp);
					}					
					//System.out.println();
					double[] normalizedForwardDensityWeights =forwardDensityWeights.clone();
					NumUtils.expNormalize(normalizedForwardDensityWeights);
					int sampledIndx = SampleUtils.sampleMultinomial(rand, normalizedForwardDensityWeights);

					logWeights[0] = forwardDensityWeights[sampledIndx];
					//					System.out.println("sampledIndx "+sampledIndx);
					//					System.out.println("weight 0: "+logWeights[0]+" normalized: "+normalizedForwardDensityWeights[sampledIndx]);
					PartialCoalescentState4BackForwardKernel2 newAncestor = (PartialCoalescentState4BackForwardKernel2) samples
							.get(sampledIndx);					
					//		double param0 = 0.1
					//				/ BackForwardKernel2.nChoose2(conditionedState.parentState().getCurrentState().nRoots());
					//		final double delta0 = newAncestor.getDelta();
					//			double logExpDensityDeltaOld = Sampling.exponentialLogDensity(param0, delta0);
					//							double refWeight =conditionedState.getCurrentState().logLikelihoodRatio()
					//							+ conditionedState.getMidState().logLikelihoodRatio()
					//							- newAncestor.getCurrentState().logLikelihoodRatio() - logExpDensityDeltaOld;
					//							System.out.println("ref weight: "+refWeight);							
					samples.set(0, samples.get(sampledIndx));							
					conditionedState.setParent(newAncestor);														
				}
				samples = resample(samples, normalizedWeights, rand);
				for(int i=0;i<normalizedWeights.length;i++) normalizedWeights[i]=1.0/N;
			}
			if (verbose)
				LogInfo.end_track();
			if (t == T - 1 && schedule != null) {
				if (resampleLastRound) {
					// this might be useful when processing a lot of particles
					// is expensive
					Pair<List<S>, double[]> resampled = resampleAndPack(samples, normalizedWeights, rand);
					samples = resampled.getFirst();
					normalizedWeights = resampled.getSecond();
				}
				if (verbose)
					LogInfo.track("Processing particles");
				for (int n = 0; n < normalizedWeights.length; n++) {
					if (verbose)
						LogInfo.logs("Particle " + (n + 1) + "/" + normalizedWeights.length);
					if (samples.get(n) != null) {
						processor.process(samples.get(n), normalizedWeights[n]);
						samples.set(n, null);
					}
				}
				if (verbose)
					LogInfo.end_track();
			}
			if (schedule != null)
				schedule.monitor(new ProcessScheduleContext(t, t == T - 1, ResampleStatus.NA));
		}
		setUnconditional();
	}

	// public static int _randcoutn = 0;

	private boolean hasNulls(List<S> samples) {
		for (S item : samples)
			if (item == null)
				return true;
		return false;
	}

	private <S> Pair<List<S>, double[]> resampleAndPack(List<S> list, double[] w, Random rand) {
		if (!NumUtils.normalize(w))
			throw new MeasureZeroException();

		if (list.size() != w.length)
			throw new RuntimeException();

		Counter<Integer> packed = Sampling.efficientMultinomialSampling(rand, w, w.length);

		double[] resultWeight = new double[packed.size()];
		List<S> resultList = new ArrayList<S>(packed.size());

		int i = 0;
		for (int itemIdx : packed.keySet()) {
			S item = list.get(itemIdx);
			resultList.add(item);
			resultWeight[i++] = packed.getCount(itemIdx);
		}

		return Pair.makePair(resultList, resultWeight);


	}

	// TODO: can be more efficient, nlog(n) instead of n^2
	private List<S> resample(final List<S> list, final double[] w, Random rand) {
		if (!NumUtils.normalize(w))
			throw new MeasureZeroException();

		if (list.size() != w.length)
			throw new RuntimeException();

		Counter<Integer> packed = Sampling.efficientMultinomialSampling(rand, w, w.length);
		final List<S> result = new ArrayList<S>(list.size());
		// System.out.print("resampled(");
		for (int itemIdx : packed.keySet()) {
			S item = list.get(itemIdx);
			// System.out.print( "" + ((PartialCoalescentState) item).id + "x" +
			// packed.getCount(itemIdx) + " ");
			for (int cur = 0; cur < packed.getCount(itemIdx); cur++)
				result.add(item);
		}
		// System.out.print(")\n");
		return result;

		// List<Integer> newAncestors = new ArrayList<Integer>(N);
		// final List<S> result = new ArrayList<S>(list.size());
		// Set<Integer> indices = new HashSet<Integer>();
		// for (int n = 0; n < list.size(); n++)
		// {
		// int index;
		// if (n == 0 && isConditional())
		// index = 0;
		// else
		// {
		// try {
		// index = SampleUtils.sampleMultinomial(rand, w); }
		// catch (RuntimeException re)
		// {
		// throw new MeasureZeroException();
		// }
		// }
		// newAncestors.add(ancestors.get(index));
		// indices.add(index);
		// result.add(list.get(index));
		// }
		// if (verbose) LogInfo.logs("EffectiveSamplingSize="+indices.size());
		// ancestors = newAncestors;
		// return result;
	}

	/**
	 * Weight is NOT in log scale, in contrast to ParticleKernel.next()
	 * 
	 * @author bouchard
	 *
	 * @param <S>
	 */
	public static interface ParticleProcessor<S> {
		public void process(S state, double weight);
	}

	public static class StoreProcessor<S> implements ParticleProcessor<S> {
		public List<S> particles = CollUtils.list();
		public List<Double> ws = CollUtils.list();

		@Override
		public void process(S state, double weight) {
			particles.add(state);
			ws.add(weight);
		}

		public S sample(Random rand) {
			final int idx = Sampling.sample(rand, ws);
			return particles.get(idx);
		}

		public S argmax() {
			int argmax = -1;
			double value = Double.NEGATIVE_INFINITY;
			for (int i = 0; i < ws.size(); i++)
				if (ws.get(i) > value) {
					value = ws.get(i);
					argmax = i;
				}
			return particles.get(argmax);
		}
	}

	public static class DoNothingProcessor<S> implements ParticleProcessor<S> {
		public void process(S state, double w) {
		}
	}

	public static class ForkedProcessor<S> implements ParticleProcessor<S> {
		public List<ParticleProcessor<S>> processors = new ArrayList<ParticleProcessor<S>>();

		@SuppressWarnings("unchecked")
		public ForkedProcessor(ParticleProcessor<S>... items) {
			processors = new ArrayList(Arrays.asList(items));
		}

		public ForkedProcessor(Collection<ParticleProcessor<S>> items) {
			this.processors = CollUtils.list(items);
		}

		public void process(S state, double weight) {
			for (ParticleProcessor<S> processor : processors)
				processor.process(state, weight);
		}
	}

	/**
	 * Given a function f and a loss l, returns min_x E l(f(x), f(X)), where the
	 * expectation is approximated using the samples from the SMC, and also min_x is
	 * restricted to the set of samples returned by SMC
	 * 
	 * @author bouchard
	 *
	 * @param <D>
	 * @param <I>
	 *            I should have .equals() implemented e.g. D =
	 *            PartialCoalescentState, I = Set <Set < Languages>> i.e. clade
	 *            representation of the forest
	 */
	public static class ParticleMapperProcessor<D, I> implements ParticleProcessor<D> {
		private final Fct<D, I> prj;
		private final Counter<I> counter = new Counter<I>();

		public static <S> ParticleMapperProcessor<S, S> saveParticlesProcessor() {
			return new ParticleMapperProcessor<S, S>(new Id<S>());
		}

		public static ParticleMapperProcessor<PartialCoalescentState, PartialCoalescentState> saveCoalescentParticlesProcessor() {
			return saveParticlesProcessor();
		}

		public ParticleMapperProcessor(Fct<D, I> prj) {
			this.prj = prj;
		}

		public void process(D state, double weight) {
			counter.incrementCount(prj.evalAt(state), weight);
		}

		public I centroid(LossFct<I> loss) {
			return new BayesRiskMinimizer<I>(loss).findMin(counter);
		}

		public String printweights() {
			String s = "";
			for (I key : counter) {
				double value = counter.getCount(key);
				s += value + ",";
			}
			return s;
		}

		public I map() {
			return counter.argMax();
		}

		public I sample(Random rand) {
			double[] probs = new double[counter.size()];
			ArrayList<I> states = new ArrayList<I>();
			int i = 0;
			for (I key : counter.keySet()) {
				states.add(key);
				probs[i++] = counter.getCount(key);
			}
			int j = Multinomial.sample(rand, probs);
			return states.get(j);
		}

		public Counter<I> getCounter() {
			return counter;
		}

	}

	public static class PCSHash implements ParticleProcessor<PartialCoalescentState> {
		private Hasher hasher = new Hasher();

		public void process(PartialCoalescentState state, double weight) {
			hasher.add(weight).add(state.logLikelihood()).add(state.topHeight())
			.add(state.getUnlabeledArbre().deepToLispString());
		}

		public int getHash() {
			return hasher.hashCode();
		}
	}

	public static class MAPDecoder<D> implements ParticleProcessor<D> {
		private D argmax = null;
		private double max = Double.NEGATIVE_INFINITY;

		public void process(D state, double weight) {
			if (weight > max) {
				argmax = state;
				max = weight;
			}
		}

		public D map() {
			return argmax;
		}
	}

}
