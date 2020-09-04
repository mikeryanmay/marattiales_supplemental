#include <Rcpp.h>
using namespace Rcpp;

// code for simulating taxon data under
// the fossilized birth-death model

// [[Rcpp::export]]
List simulateEpochsCPP(double age,
                       NumericVector epoch_times,
                       NumericVector lambda,
                       NumericVector mu,
                       NumericVector phi,
                       double rho) {

    // initialize some variables
    size_t num_epochs = epoch_times.size();
    int current_num_taxa = 1;
    NumericVector num_taxa_per_epoch(num_epochs + 1); // (the first value is for the present)

    // initialize LTT parameters
    std::vector<int>    num_taxa_per_time(0);
    std::vector<double> event_times(0);
    
    num_taxa_per_time.push_back(current_num_taxa);
    event_times.push_back(age);
    
    // initialize the current state
    double current_time = age;
    double next_time;

    // get the current epoch
    int current_epoch = 0;
    for(size_t i = 0; i < num_epochs; ++i) {
        if ( epoch_times[i] > age ) {
            break;
        } else {
            current_epoch++;
        }
    }

    // compute times of the next epoch
    int next_epoch         = current_epoch == 0 ? 0 : current_epoch - 1;
    double next_epoch_time = epoch_times[next_epoch];

    // Rcout << current_epoch << " -- " << current_time << std::endl;
    // Rcout << next_epoch << " -- " << next_epoch_time << std::endl;

    // get the current rates
    double current_lambda = lambda[current_epoch];
    double current_mu     = mu[current_epoch];
    double current_phi    = phi[current_epoch];
    double total_rate     = current_lambda + current_mu + current_phi;

    bool end = false;
    while (end == false) {

        // draw the next time
        next_time = current_time - rexp(1, double(current_num_taxa) * total_rate)[0];

        if ( next_time < next_epoch_time ) {

            if ( next_epoch < 0 ) {

                // terminate if the next epoch is the present
                end = true;

            } else {

                // we moved into the next epoch
                current_time    = next_epoch_time;
                current_epoch   = next_epoch;
                next_epoch      = current_epoch - 1;

                // compute the rates
                current_lambda = lambda[current_epoch];
                current_mu     = mu[current_epoch];
                current_phi    = phi[current_epoch];
                total_rate     = current_lambda + current_mu + current_phi;

                // make sure the next time is up-to-date
                if ( next_epoch < 0 ) {
                    next_epoch_time = 0;
                } else {
                    next_epoch_time = epoch_times[next_epoch];
                }

            }


        } else {

            // increment time
            current_time = next_time;

            // simulate an event
            double u = runif(1, 0, total_rate)[0];

            if ( u < current_lambda ) { // speciation event
                current_num_taxa++;
                num_taxa_per_time.push_back(current_num_taxa);
                event_times.push_back(current_time);
            } else if ( u < current_lambda + current_mu ) { // extinction event
                current_num_taxa--;
                num_taxa_per_time.push_back(current_num_taxa);
                event_times.push_back(current_time);
            } else {
                num_taxa_per_epoch[current_epoch+1]++;
            }

        }

        // Rcout << current_time << " -- " << current_num_taxa << std::endl;

        // check for extinction
        if (current_num_taxa == 0) {
            end = true;
        }

        // if ( current_num_taxa * rho > 1000 ) {
        //     current_num_taxa = 0;
        //     end = true;
        // }

    }

    // increment the events
    current_time = current_time < 0 ? 0 : current_time;
    num_taxa_per_time.push_back(current_num_taxa);
    event_times.push_back(current_time);
    
    if ( current_num_taxa == 0 ) {
        return Rcpp::List::create(Rcpp::Named("num_taxa_per_epoch") = num_taxa_per_epoch,
                                  Rcpp::Named("num_taxa_per_time")  = num_taxa_per_time,
                                  Rcpp::Named("event_times")        = event_times);
    }

    // simulate taxon sampling
    double new_num_taxa = rbinom(1, current_num_taxa, rho)[0];

    // append extant taxa
    if ( current_num_taxa > 0 ) {
        num_taxa_per_epoch[0] = new_num_taxa;
    }

    return Rcpp::List::create(Rcpp::Named("num_taxa_per_epoch") = num_taxa_per_epoch,
                              Rcpp::Named("num_taxa_per_time")  = num_taxa_per_time,
                              Rcpp::Named("event_times")        = event_times);
    
    // return num_taxa_per_epoch;

}
















