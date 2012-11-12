#ifndef BESTPROBABILITY_H
#define BESTPROBABILITY_H

class BestProbability
{
  template<class Chart>
  static void calculate_maxrule_probabilities(Chart & chart);
};




template<class Chart>
void
BestProbability::calculate_maxrule_probabilities (Chart & chart)
{
    chart->opencells_apply_bottom_up(
      [](typename Chart::Cell & cell){ cell.calculate_maxrule_probabilities(); }
  );
}


#endif