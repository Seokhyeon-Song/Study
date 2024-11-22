#include <iostream>
#include <fstream>
#include <cmath>
#include <complex>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>

gsl_complex Pi(double p, double Ms, double M) {
	double zA = 0, zB = 0;
	double sqfactor = 0;
	switch (gsl_fcmp(2 * Ms, M, 0.000001)) {
	case 0:	// 2 * Ms == M
		zA = 2 / gsl_pow_2(M);
		zB = -2;
		break;
	case 1:	// 2 * Ms > M
		sqfactor = std::sqrt(4 * gsl_pow_2(Ms) - gsl_pow_2(M));
		zA = (1 / gsl_pow_2(M)) - (4 * gsl_pow_2(Ms) / (sqfactor * gsl_pow_3(M))) * std::atan(M / sqfactor);
		zB = (12 * gsl_pow_2(Ms) - 2 * gsl_pow_2(M)) * std::atan(M / sqfactor) / (sqfactor * M) - 1;
		break;
	case -1: // 2 * Ms < M
		sqfactor = std::sqrt(gsl_pow_2(M) - 4 * gsl_pow_2(Ms));
		zA = (1 / gsl_pow_2(M)) + (4 * gsl_pow_2(Ms) / (sqfactor * gsl_pow_3(M))) * std::atanh(sqfactor / M);
		zB = (2 * gsl_pow_2(M) - 12 * gsl_pow_2(Ms)) * std::atanh(sqfactor / M) / (sqfactor * M) - 1;
		break;
	}

	gsl_complex pibare;
	switch (gsl_fcmp(2 * Ms, p, 0.000001)) {
	case 0:
		GSL_SET_COMPLEX(&pibare, 0.0, 0.0);
		break;
	default:
		gsl_complex sqfactor2;
		GSL_SET_COMPLEX(&sqfactor2, 4 * gsl_pow_2(Ms) - gsl_pow_2(p), 0.0);
		gsl_complex atanfactor = gsl_complex_mul_real(gsl_complex_inverse(gsl_complex_sqrt(sqfactor2)), p);
		pibare = gsl_complex_mul_real(gsl_complex_div(gsl_complex_arctan(atanfactor), atanfactor), -2);
		break;
	}

	return gsl_complex_add_real(pibare, zA * gsl_pow_2(p) + zB);
}

struct null_hypo_params {
	double A;
	double M;
	double Gamma;
};

struct alt_hypo_params {
	double A;
	double Msmall;
	double M;
	double Gamma;
	double C;
};

struct null_alt_params {
	null_hypo_params npars;
	alt_hypo_params apars;
};

double null_hypo_dist(double p, null_hypo_params pars)
{
	return pars.A * gsl_pow_2(p) / (gsl_pow_2(gsl_pow_2(p) - gsl_pow_2(pars.M)) + gsl_pow_2(pars.Gamma));
}

double alt_hypo_dist(double p, alt_hypo_params pars)
{
	gsl_complex tree, iGamma;
	GSL_SET_COMPLEX(&tree, gsl_pow_2(p) - gsl_pow_2(pars.M), 0.0);
	GSL_SET_COMPLEX(&iGamma, 0.0, pars.Gamma);
	return pars.A * gsl_pow_2(p) / gsl_complex_abs2(gsl_complex_add(gsl_complex_add(tree, iGamma), gsl_complex_mul_real(Pi(p, pars.Msmall, pars.M), pars.C)));
}

double null_hypo_loglike(double p, void* theta) {
	null_alt_params& pars = *reinterpret_cast<null_alt_params*>(theta);
	double null_prob = null_hypo_dist(p, pars.npars);
	double alt_prob = alt_hypo_dist(p, pars.apars);
	return alt_prob * std::log(null_prob) + (1 - alt_prob) * std::log(1 - null_prob);
}

double alt_hypo_loglike(double p, void* theta) {
	alt_hypo_params& pars = *reinterpret_cast<alt_hypo_params*>(theta);
	double prob = alt_hypo_dist(p, pars);
	return prob * std::log(prob) + (1 - prob) * std::log(1 - prob);
}

double graddesc_A(double p, void* theta) {
	null_alt_params& pars = *reinterpret_cast<null_alt_params*>(theta);
	double null_prob = null_hypo_dist(p, pars.npars);
	double alt_prob = alt_hypo_dist(p, pars.apars);
	double partial_A = gsl_pow_2(p) / (gsl_pow_2(gsl_pow_2(p) - gsl_pow_2(pars.npars.M)) + gsl_pow_2(pars.npars.Gamma));
	return ((alt_prob / null_prob) - ((1 - alt_prob) / (1 - null_prob))) * partial_A;
}

double graddesc_M(double p, void* theta) {
	null_alt_params& pars = *reinterpret_cast<null_alt_params*>(theta);
	double null_prob = null_hypo_dist(p, pars.npars);
	double alt_prob = alt_hypo_dist(p, pars.apars);
	double partial_M = 4 * pars.npars.M * (gsl_pow_2(p) - gsl_pow_2(pars.npars.M)) * pars.npars.A * gsl_pow_2(p) / gsl_pow_2(gsl_pow_2(gsl_pow_2(p) - gsl_pow_2(pars.npars.M)) + gsl_pow_2(pars.npars.Gamma));
	return ((alt_prob / null_prob) - ((1 - alt_prob) / (1 - null_prob))) * partial_M;
}

double graddesc_Gamma(double p, void* theta) {
	null_alt_params& pars = *reinterpret_cast<null_alt_params*>(theta);
	double null_prob = null_hypo_dist(p, pars.npars);
	double alt_prob = alt_hypo_dist(p, pars.apars);
	double partial_Gamma = - 2 * pars.npars.Gamma * pars.npars.A * gsl_pow_2(p) / gsl_pow_2(gsl_pow_2(gsl_pow_2(p) - gsl_pow_2(pars.npars.M)) + gsl_pow_2(pars.npars.Gamma));
	return ((alt_prob / null_prob) - ((1 - alt_prob) / (1 - null_prob))) * partial_Gamma;
}

double alt_hypo_signalprob(double p, void* theta) {
	alt_hypo_params& pars = *reinterpret_cast<alt_hypo_params*>(theta);
	return alt_hypo_dist(p, pars);
}

struct fitresult {
	double A;
	double M;
	double Gamma;
	double NS;
};

fitresult bestfit(alt_hypo_params apars, gsl_integration_workspace* workspace) {
	// Integration Setup
	double totalwidth = apars.Gamma + apars.C * GSL_IMAG(Pi(apars.M, apars.Msmall, apars.M));
	double result_A, result_M, result_Gamma, result_n, prev_result_n, result_a, error;
	const double xlow = apars.M - 2 * totalwidth;
	const double xhigh = apars.M + 2 * totalwidth;

	// Initial Guess
	null_hypo_params npars;
	npars.A = apars.A;
	npars.M = apars.M;
	npars.Gamma = totalwidth;

	// Gradient.. Ascent?
	double learning_rate = 1;
	null_alt_params napars = { npars, apars };
	null_alt_params prev_napars = napars;

	double add_A = 0.5 * npars.A, add_M = 0.5 * npars.Gamma, add_Gamma = 0.5 * npars.Gamma;

	gsl_function gdesc_A, gdesc_M, gdesc_Gamma, nloglike, aloglike;
	gdesc_A.function = &graddesc_A;
	gdesc_M.function = &graddesc_M;
	gdesc_Gamma.function = &graddesc_Gamma;
	nloglike.function = &null_hypo_loglike;
	aloglike.function = &alt_hypo_loglike;
	gdesc_A.params = reinterpret_cast<void*>(&napars);
	gdesc_M.params = reinterpret_cast<void*>(&napars);
	gdesc_Gamma.params = reinterpret_cast<void*>(&napars);
	nloglike.params = reinterpret_cast<void*>(&napars);
	gsl_integration_qag(&nloglike, xlow, xhigh, 1e-14, 1e-6, 5000, 6, workspace, &prev_result_n, &error);

	for (int i = 0; i < 200; i++) {
		gsl_integration_qag(&gdesc_A, xlow, xhigh, 1e-10, 1e-5, 5000, 6, workspace, &result_A, &error);
		gsl_integration_qag(&gdesc_M, xlow, xhigh, 1e-10, 1e-5, 5000, 6, workspace, &result_M, &error);
		gsl_integration_qag(&gdesc_Gamma, xlow, xhigh, 1e-10, 1e-5, 5000, 6, workspace, &result_Gamma, &error);

		while (true) {
			napars = prev_napars;
			napars.npars.A += add_A * GSL_SIGN(result_A);
			gsl_integration_qag(&nloglike, xlow, xhigh, 1e-14, 1e-6, 5000, 6, workspace, &result_n, &error);
			if (result_n > (1 + 1e-6) * prev_result_n)
				break;
			add_A *= 0.5;
		}
		while (true) {
			napars = prev_napars;
			napars.npars.M += add_M * GSL_SIGN(result_M);
			gsl_integration_qag(&nloglike, xlow, xhigh, 1e-14, 1e-6, 5000, 6, workspace, &result_n, &error);
			if (result_n > (1 + 1e-6) * prev_result_n)
				break;
			add_M *= 0.5;
		}
		while (true) {
			napars = prev_napars;
			napars.npars.Gamma += add_Gamma * GSL_SIGN(result_Gamma);
			gsl_integration_qag(&nloglike, xlow, xhigh, 1e-14, 1e-6, 5000, 6, workspace, &result_n, &error);
			if (result_n > (1 + 1e-6) * prev_result_n)
				break;
			add_Gamma *= 0.5;
		}
		while (true) {
			napars = prev_napars;
			napars.npars.A += add_A * GSL_SIGN(result_A);
			napars.npars.M += add_M * GSL_SIGN(result_M);
			napars.npars.Gamma += add_Gamma * GSL_SIGN(result_Gamma);
			gsl_integration_qag(&nloglike, xlow, xhigh, 1e-14, 1e-6, 5000, 6, workspace, &result_n, &error);
			if (result_n > (1 + 1e-6) * prev_result_n) {
				prev_result_n = result_n;
				break;
			}
			add_A *= 0.5;
			add_M *= 0.5;
			add_Gamma *= 0.5;
		}

		if (add_A / napars.npars.A < 1e-6 && add_M / napars.npars.Gamma < 1e-6 && add_Gamma / napars.npars.Gamma < 1e-6) {
			break;
		}
		prev_napars = napars;
	}

	aloglike.params = reinterpret_cast<void*>(&apars);
	gsl_integration_qag(&aloglike, xlow, xhigh, 1e-14, 1e-6, 5000, 6, workspace, &result_a, &error);

	gsl_function aprob;
	aprob.function = &alt_hypo_signalprob;
	aprob.params = reinterpret_cast<void*>(&apars);
	double asignalprob;
	gsl_integration_qag(&aprob, xlow, xhigh, 1e-14, 1e-6, 5000, 6, workspace, &asignalprob, &error);

	double NS = 2 * totalwidth * asignalprob * 23.66 / (result_a - result_n);

	return { napars.npars.A, napars.npars.M, napars.npars.Gamma, NS };
}

int main() {
	gsl_set_error_handler_off();
	std::ofstream writeFile;
	writeFile.open("fitresult.txt");

	gsl_integration_workspace* workspace = gsl_integration_workspace_alloc(10000);
	// True Parameters
	alt_hypo_params apars;
	fitresult res;
	apars.M = 1;
	apars.C = 0.04;

	/*
	//Debug
	apars.Gamma = 0.103;
	apars.Msmall = 0.41966;
	double totalwidth = apars.Gamma + apars.C * GSL_IMAG(Pi(apars.M, apars.Msmall, apars.M));
	apars.A = 1e-2 * gsl_pow_2(totalwidth) / gsl_pow_2(apars.M);
	res = bestfit(apars, workspace);
	std::cout << res.NS << std::endl;
	*/
	
	int domain_resolution = 100;
	writeFile << "alpha/4Pi: " << apars.C << std::endl;
	writeFile << "Mother particle mass is fixed to M=1GeV" << std::endl;
	writeFile << "vertical: Gamma = " << 0.4 / domain_resolution << "M to 0.4M" << std::endl;
	writeFile << "horizontal: m = 0.5M-Gamma to 0.5M+Gamma" << std::endl;

	for (int gamma_div = 1; gamma_div <= domain_resolution; gamma_div++) {
		for (int m_div = 0; m_div <= domain_resolution; m_div++) {
			apars.Gamma = 0.1 * apars.M + 0.3 * apars.M * gamma_div / domain_resolution;
			apars.Msmall = 0.5 * apars.M - apars.Gamma + 2.0 * apars.Gamma * m_div / domain_resolution;
			double totalwidth = apars.Gamma + apars.C * GSL_IMAG(Pi(apars.M, apars.Msmall, apars.M));
			apars.A = 1e-3 * gsl_pow_2(totalwidth) / gsl_pow_2(apars.M);

			res = bestfit(apars, workspace);
			writeFile << (int)res.NS << " ";
			std::cout << m_div << std::endl;
		}
		writeFile << std::endl;
		std::cout << 100.0 * gamma_div / domain_resolution << "% done..." << std::endl;
	}
	
	
	writeFile.close();
	return 0;
}