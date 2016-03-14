#ifndef CTCLMI_H
#define	CTCLMI_H

#include <ibex/ibex.h>
#include <sdpa_call.h>

class CtcLMI : public ibex::Ctc {
public:

    CtcLMI(ibex::MatrixArray &matrices) : Ctc(-1), matrices(matrices) {
    };
    void contract(ibex::IntervalVector& box);

private:
    SDPA solver;
    void mkSolver(ibex::MatrixArray &matrices, const ibex::IntervalVector &box);
    ibex::MatrixArray &matrices;
};

#endif	/* CTCLMI_H */

