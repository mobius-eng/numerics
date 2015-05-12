/*
 * Copyright (c) 2015 Alexey Cherkaev.
 * This file is licensed under Lesser General Public License v.3
 * The text of the license can be found at http://www.gnu.org/licenses/lgpl-3.0.txt
 */

package numerics;

import clojure.lang.AFn;
import clojure.lang.PersistentVector;


public final class Polint {

    private static int findClosestIndex(double[] v, double x) {
        int len = v.length;
        int n;
        double diffn;
        int i = 0;
        double diffi = Math.abs(x - v[i]);
        int j = len-1;
        double diffj = Math.abs(x - v[j]);
        for(;;) {
            if (Math.abs(i-j) <= 1) {
                break;
            }
            else {
                n = (i+j) / 2;
                diffn = Math.abs(x - v[n]);
                if (diffn == 0.0) {
                    diffi = 0.0;
                    diffj = 0.0;
                    i = n;
                    j = n;
                    break;
                }
                if (Math.abs(diffn - diffi) > Math.abs(diffn - diffj)) {
                    i = n;
                    diffi = diffn;
                }
                else {
                    j = n;
                    diffj = diffn;
                }
            }
        }
        return ((diffi < diffj) ? i : j);
    }

    public static AFn neville (final double[] xa, final double[] ya) {
        final int n = xa.length;
        return new AFn() {
            @Override public Object invoke(Object arg) {
                double x = (double) arg;
                double coeff;
                double xi_x;
                double xim_x;
                double cd_diff;
                double y;
                double dy = 0.0;
                double[] c = ya.clone();
                double[] d = ya.clone();
                int ind = findClosestIndex(xa, x);
                y = ya[ind];
                for (int m = 1; m < n; m++) {
                    for (int i = 0; i < n - m; i++) {
                        xi_x = xa[i] - x;
                        xim_x = xa[i + m] - x;
                        cd_diff = c[i + 1] - d[i];
                        coeff = cd_diff / (xa[i] - xa[i + m]);
                        d[i] = xim_x * coeff;
                        c[i] = xi_x * coeff;
                    }
                    if (2 * ind < n - m) {
                        dy = c[ind];
                    }
                    else {
                        dy = d [(ind-1)];
                        ind--;
                    }
                    y += dy;
                }
                return PersistentVector.create(y, dy);
            }
        };
    }
}
