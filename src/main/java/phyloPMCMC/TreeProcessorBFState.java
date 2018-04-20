package phyloPMCMC;


import smc.PartialCoalescentState4BackForwardKernel;
import ev.poi.processors.TreeDistancesProcessor;

public class TreeProcessorBFState extends TreeDistancesProcessor
  
{ 
  @Override
  public void process(Object state, double weight)
  {
    if (state instanceof PartialCoalescentState4BackForwardKernel)
        super.process(((PartialCoalescentState4BackForwardKernel) state).getCurrentState(), weight);    
 
  }
 
}
