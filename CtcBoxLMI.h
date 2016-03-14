#ifndef CTCBOXLMI_H
#define	CTCBOXLMI_H

#include <ibex/ibex.h>
#include <sdpa_call.h>

class CtcBoxLMI : public ibex::Ctc {
public:

    CtcBoxLMI(const ibex::MatrixArray &mA) : Ctc(-1), matrices(mA) {
    };

    void contract(ibex::IntervalVector& box);
private:
    SDPA solver;
    void mkSolver(const ibex::MatrixArray &matrices, ibex::IntervalVector &box);
    const ibex::MatrixArray &matrices;
};

#endif	/* CTCBOXLMI_H */

