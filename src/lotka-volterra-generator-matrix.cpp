#include <Rcpp.h>
#include <vector>
using namespace std;
using namespace Rcpp;


vector<int> obtain_state_counts_lotka_volterra(int y1, int y2, int state, int y1_min, int y1_max, int y2_min, int y2_max){

  vector<int> ys(2);

  ys[0] = (int)floor((state-1)/(1.0*(y2_max-y2_min+1))) + y1_min;
  ys[1] = state - (ys[0] - y1_min)*(y2_max - y2_min + 1) + y2_min - 1;

  return(ys);

}

int obtain_state_lotka_volterra(int y1, int y2, int y1_min, int y1_max, int y2_min, int y2_max){

  return ((y1 - y1_min)*(y2_max - y2_min + 1) + (y2 - y2_min + 1));

}

// [[Rcpp::export]]
List lotka_volterra_generator(NumericVector theta, int y1_min, int y1_max, int y2_min, int y2_max) {

  double cum_sum;
  int num_states, i, state;
  vector<int> ys(2);

  num_states = (y1_max-y1_min+1)*(y2_max-y2_min+1);

  // dynamic vector for rows, cols, values:
  vector<int> rows;
  vector<int> cols;
  vector<double> values;

  for (i = 1; i<=num_states; i++){
    cum_sum = 0;

    ys = obtain_state_counts_lotka_volterra(ys[0], ys[1], i, y1_min, y1_max, y2_min, y2_max);

    // check to see if a prey can die
    if (ys[0] < y1_max){
      state = obtain_state_lotka_volterra(ys[0]+1,ys[1],y1_min,y1_max,y2_min,y2_max);
      rows.push_back(i);
      cols.push_back(state);
      values.push_back(theta[0]*ys[0]);
      cum_sum += values.back(); // get last value, i.e. theta[0]*ys[0]
    }

    // check to see if predator can kill a prey
    if (ys[0] > y1_min && ys[1]< y2_max){
      state = obtain_state_lotka_volterra(ys[0]-1,ys[1]+1,y1_min,y1_max,y2_min,y2_max);
      rows.push_back(i);
      cols.push_back(state);
      values.push_back(theta[1]*ys[0]*ys[1]);
      cum_sum += values.back();
    }

    // check to see if a predator can die
    if (ys[1]> y2_min){
      state = obtain_state_lotka_volterra(ys[0], ys[1]-1, y1_min, y1_max, y2_min, y2_max);
      rows.push_back(i);
      cols.push_back(state);
      values.push_back(theta[2]*ys[1]);
      cum_sum += values.back();
    }

    // now make sure we get the diagonal term correct
    rows.push_back(i);
    cols.push_back(i);
    values.push_back(-cum_sum);

  }

  return List::create(Named("rows") = wrap(rows),
                      Named("cols") = wrap(cols),
                      Named("vals") = wrap(values) );

}

