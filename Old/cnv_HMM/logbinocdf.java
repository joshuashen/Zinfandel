package cnv_HMM;
import cnv_simulation.*;
import org.apache.commons.math.util.ContinuedFraction;
import org.apache.commons.math.special.Beta;
import org.apache.commons.math.MathException;

public class logbinocdf {
    static double epsilon = 1e-14;

    static double compute(int X, int N, double P)
    throws MathException {
        if (X < 0) {
            return Double.NEGATIVE_INFINITY;
        } else if (X >= N || P <= (X+2.)/(N+3))  {
            return 0;
        } else {
            return
                logRegularizedBeta(1-P, N-X, X+1);
        }
    }

    static double logRegularizedBeta(double x, final double a, final double b)
    throws MathException {
            ContinuedFraction fraction = new ContinuedFraction() {

                protected double getB(int n, double x) {
                    double ret;
                    double m;
                    if (n % 2 == 0) { // even
                        m = n / 2.0;
                        ret = (m * (b - m) * x) /
                            ((a + (2 * m) - 1) * (a + (2 * m)));
                    } else {
                        m = (n - 1.0) / 2.0;
                        ret = -((a + m) * (a + b + m) * x) /
                                ((a + (2 * m)) * (a + (2 * m) + 1.0));
                    }
                    return ret;
                }

                protected double getA(int n, double x) {
                    return 1.0;
                }
            };
            double ret1 = (a * Math.log(x)) + (b * Math.log(1.0 - x)) - Math.log(a);
	    double ret2 = Beta.logBeta(a, b, epsilon, Integer.MAX_VALUE);
            double ret3 = fraction.evaluate(x, epsilon, Integer.MAX_VALUE);
	    double ret4 = Math.log(ret3);
	    if(Double.isNaN(ret1) || Double.isNaN(ret2) || Double.isNaN(ret4))
		throw new MathException("Bad calculation for x " + x + " a " + a + " b " + b + " ret1 " + ret1 + " ret2 " + ret2 + " ret3 " + ret3);
	    return ret1 - ret2 - ret4;
    }
}
