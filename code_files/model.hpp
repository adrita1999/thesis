/*****************************************************************************
 *   kmdiff
 *   Authors: T. Lemane
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Affero General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Affero General Public License for more details.
 *
 *  You should have received a copy of the GNU Affero General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *****************************************************************************/

#pragma once

#include <algorithm>
#include <cmath>
#include <numeric>
#include <random>
#include <vector>


#include <alglibinternal.h>
#include <spdlog/spdlog.h>
#include <specialfunctions.h>

#define KMTRICKS_PUBLIC
#include <kmtricks/utils.hpp>

#include <kmdiff/range.hpp>
#include <kmdiff/kmer.hpp>
#include <kmdiff/log_factorial_table.hpp>
#include <kmdiff/correction.hpp>
#include <kmdiff/imodel.hpp>
//modified
#include <mutex>
//





namespace kmdiff {

  template <size_t MAX_C>
  class Model
  {
    using count_t = typename km::selectC<MAX_C>::type;

   public:
    Model() {}

    template <typename Iterable>
    double compute_mean(Iterable&& v)
    {
      return std::accumulate(v.begin(), v.end(), 0.0) / v.size();
    }

    template <typename Iterable>
    std::tuple<double, size_t> compute_mean_e(Iterable& v)
    {
      size_t sum = 0;
      size_t positive = 0;
      for (auto& e : v)
      {
        if (e > 0) positive++;
        sum += e;
      }
      double mean = sum / static_cast<double>(v.size());
      return std::make_tuple(mean, positive);
    }

    template <typename Iterable>
    std::tuple<double, size_t> compute_sum_e(Iterable& v)
    {
      double sum = 0;
      size_t positive = 0;
      for (auto& e : v)
      {
        if (e > 0) positive++;
        sum += e;
        // std::cout<<e<< ", ";
      }
      // std::cout<<"\n";
      return std::make_tuple(sum, positive);
    }

    template <typename Iterable>
    std::tuple<double, double> compute_mean_sd(Iterable& v)
    {
      double mean = compute_mean(v);
      double sq_sum = std::inner_product(v.begin(), v.end(), v.begin(), 0.0);
      double sd = std::sqrt(sq_sum / v.size() - mean * mean);
      return std::make_tuple(mean, sd);
    }

  };

  template <size_t MAX_C>
  class PoissonLikelihood : public IModel<MAX_C>, public Model<MAX_C>
  {
  /*
    Based on https://github.com/atifrahman/HAWK
    https://doi.org/10.7554/eLife.32920.001
    https://doi.org/10.1371/journal.pone.0245058
  */

    using count_t = typename km::selectC<MAX_C>::type;

   public:
    PoissonLikelihood(
        size_t nb_controls,
        size_t nb_cases,
        const std::vector<size_t>& total_controls,
        const std::vector<size_t>& total_cases,
        size_t preload,
        //modified
        double *pheno
        //
        )
        : m_nb_controls(nb_controls),
          m_nb_cases(nb_cases),
          m_total_kmer_controls(total_controls),
          m_total_kmer_cases(total_cases),
          m_preload(preload),
          //modified
          phenotypeVals(pheno)
          //
    {
        //modified
        noSamples=m_nb_controls+m_nb_cases;
        //std::cout<<"constructor\n";
        //check
        // for (int i = 0; i < noSamples; i++)
        // {
        //   std::cout<<phenotypeVals[i]<<"  "<<&phenotypeVals[i]<<"\n";
        // }
        //
    }
    
    void configure(const std::string& config) override {}
   private:
    double log_factorial(int k)
    {
      double res = 0;
      while (k > 1)
      {
        res += log(k);
        k--;
      }
      return res;
    }

    double poisson_prob(int k, double lambda)
    {
      if (lambda <= 0) return 0;
      if (k < 0) k = 0;
      return (-lambda + (k * log(lambda) - m_lf_table[k]));
    }

   public:

    //modified
    void computeNullLikelihood() override
    {
      double y_mean = 0;

      for (int k = 0; k < noSamples; k++)
      {
        y_mean += phenotypeVals[k];
      }

      y_mean = y_mean / noSamples;

      double e_null = 0;

      for (int k = 0; k < noSamples; k++)
      {
        e_null += (phenotypeVals[k] - y_mean) * (phenotypeVals[k] - y_mean);
      }

      e_null = e_null / noSamples;
      // e_null = sqrt(e_null / noSamples); // variance of phenotype value w.r.t mean

      likelihoodNull = 0;
      double logPI=0.9189385;

      for (int k = 0; k < noSamples; k++) 
      {
        likelihoodNull+= (-log(e_null) / 2.0 - logPI - (phenotypeVals[k] - y_mean) * (phenotypeVals[k] - y_mean) / (2 * e_null));
      }
      //std::cout<<"likelihood null "<<likelihoodNull<<"\n";
    }



    std::vector<double> regress(int noSamples, double *x, double *y)
    {
      std::vector<double> result;
      double a, b;
      double x_mean = 0, y_mean = 0;
      double s_x = 0, s_y = 0, s_xy = 0;

      for (int i = 0; i < noSamples; i++)
      {
        x_mean += x[i];
        y_mean += y[i];
      }
      x_mean = x_mean / noSamples;
      y_mean = y_mean / noSamples;

      for (int i = 0; i < noSamples; i++)
      {
        s_x += (x[i] - x_mean) * (x[i] - x_mean);
        s_y += (y[i] - y_mean) * (y[i] - y_mean);
        s_xy += (x[i] - x_mean) * (y[i] - y_mean);
      }

      b = s_xy / s_x;			 // slope
      a = y_mean - b * x_mean; // y-intercept
      // regression using covariance and correlation between ith kmer count and phenotype value
      result.push_back(a);
      result.push_back(b);
      return result;
    }


    //

    std::tuple<double, Significance, double, double>
    process(const Range<count_t>& controls, const Range<count_t>& cases) override
    {
      //std::cout<<"model.hpp -> poisson er process call hocche\n";
      // std::cout<<"case "<<m_total_kmer_cases.size()<<"\n";
      // std::cout<<"control "<<m_total_kmer_controls.size()<<"\n";

      auto [mean_control, positive_controls] = this->compute_sum_e(controls);
      auto [mean_case, positive_cases] = this->compute_sum_e(cases);

      double mean = (mean_control + mean_case) / static_cast<double>(m_sum_controls + m_sum_cases);

      // prev code

      // double null_hypothesis = 0;
      // double alt_hypothesis = 0;

      // alt_hypothesis += poisson_prob(mean_control, mean_control);
      // alt_hypothesis += poisson_prob(mean_case, mean_case);

      // null_hypothesis += poisson_prob(mean_control, mean * m_sum_controls);
      // null_hypothesis += poisson_prob(mean_case, mean * m_sum_cases);

      // double likelihood_ratio = alt_hypothesis - null_hypothesis;

      // if (likelihood_ratio < 0) likelihood_ratio = 0;
      // double p_value = alglib::chisquarecdistribution(1, 2 * likelihood_ratio);

      //

      //modified

      double *x=new double[noSamples];


      // for(int i=0;i<m_nb_controls;i++)
      // {
      //   std::cout<<"control "<<i<<" "<<controls[i]<<"\n";
      // }

      // for (int i=0;i<noSamples;i++)
      // {
      //   if(i<m_nb_controls)
      //   {
      //       x[i]=controls[i]/m_total_kmer_controls[i];
      //   }
      //   else
      //   {
      //       x[i]=cases[i-m_nb_controls]/m_total_kmer_controls[i-m_nb_controls];
      //   }
      // }

      int iter=0;
      for (auto& e : controls)
      {
          // mtx_print.lock();
          // std::cout<<"e "<<e<<"\n";
          // mtx_print.unlock();
          // mtx_print.lock();
          // std::cout<<"tot "<<m_total_kmer_controls[iter]<<"\n";
          // mtx_print.unlock();
          x[iter]=(e*1.0)/m_total_kmer_controls[iter];
          iter++;
      }
      for (auto& e : cases)
      {
          x[iter]=(e*1.0)/m_total_kmer_cases[iter-m_nb_controls];
          iter++;
      }


      // for(int i=0;i<noSamples;i++)
      // {
      //   mtx_print.lock();
      //   std::cout<<"x "<<i<<" "<<x[i]<<"\n";
      //   mtx_print.unlock();
      // }
      std::vector<double> result = regress(noSamples, x, phenotypeVals);
      // std::cout<<"result0 "<<result[0]<<"\n";
      // std::cout<<"result1 "<<result[1]<<"\n";


      double e_alt, y_p;
			e_alt = 0;

			for (int k = 0; k < noSamples; k++)
			{
				y_p = result[0] + result[1] * x[k];
				e_alt += (phenotypeVals[k] - y_p) * (phenotypeVals[k] - y_p);
			}

			e_alt = e_alt / noSamples;
      //std::cout<<"e_alt"<<e_alt<<"\n";
			// e_alt = sqrt(e_alt / noSamples); #variance needed, not stdev

			double likelihoodAlt = 0;
      double logPI=0.9189385;

			likelihoodAlt = (-log(e_alt) / 2.0 * noSamples - logPI * noSamples - noSamples / 2.0);
      // std::cout<<"lik null "<<likelihoodNull<<"\n"; 
      // std::cout<<"lik alt "<<likelihoodAlt<<"\n";

			double likelihoodRatio = likelihoodAlt - likelihoodNull;

			if (likelihoodRatio < 0)
			{
				likelihoodRatio = 0;
			}


      mtx_print.lock();
      std::cout<<"likelihood ratio: "<<likelihoodRatio<<"\n";
      mtx_print.unlock();

			double p_value = alglib::chisquarecdistribution(1, 2 * likelihoodRatio);

      //

      Significance sign;

      mean_control = mean_control * m_sum_cases / m_sum_controls;

      if (mean_control < mean_case)
        sign = Significance::CASE;
      else if (mean_control > mean_case)
        sign = Significance::CONTROL;
      else
        sign = Significance::NO;

      return std::make_tuple(p_value, sign, mean_control, mean_case);
    }

   private:
    size_t m_nb_controls{0};
    size_t m_nb_cases{0};

    const std::vector<size_t>& m_total_kmer_controls;
    const std::vector<size_t>& m_total_kmer_cases;

    size_t m_sum_controls{
        std::accumulate(m_total_kmer_controls.begin(), m_total_kmer_controls.end(), 0ULL)};

    size_t m_sum_cases{std::accumulate(m_total_kmer_cases.begin(), m_total_kmer_cases.end(), 0ULL)};

    size_t m_preload;
    LogFactorialTable m_lf_table{m_preload};

    //modified
    double *phenotypeVals;
    int noSamples;
    double likelihoodNull;
    std::mutex mtx_print;
    //

  };

} // end of namespace kmdiff

