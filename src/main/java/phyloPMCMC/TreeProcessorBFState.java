package phyloPMCMC;


import smc.PartialCoalescentState4BackForwardKernel2;
import ev.poi.processors.TreeDistancesProcessor;

public class TreeProcessorBFState extends TreeDistancesProcessor
  
{ 
  @Override
  public void process(Object state, double weight)
  {
    if (state instanceof PartialCoalescentState4BackForwardKernel2)
        super.process(((PartialCoalescentState4BackForwardKernel2) state).getCurrentState(), weight);    
 
  }
 
}
