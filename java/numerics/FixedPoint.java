/*
 * Copyright (c) 2015 Alexey Cherkaev.
 * This file is licensed under Lesser General Public License v.3
 * The text of the license can be found at http://www.gnu.org/licenses/lgpl-3.0.txt
 */

package numerics;

import numerics.criteria.Criteria;
import clojure.lang.IFn;
import clojure.lang.Keyword;


/**
 * Created by alex on 10/05/15.
 * Doesn't have a point: pure Clojure fixed-point is just as fast
 */
public class FixedPoint {
    public static Criteria solve(IFn checkStatus, IFn f, Object x0, Keyword cont, Keyword failed, Keyword finished) {
        IFn check = (IFn)checkStatus.invoke(x0);
        Object x = x0;
        Criteria sx = null;
        for(;;) {
            x = f.invoke(x);
            sx = (Criteria)check.invoke(x);
            if (sx.status == failed || sx.status == finished) {
                break;
            }
        }
        return sx;
    }
}
