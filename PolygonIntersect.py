#import java.awt.geom.Point2D;

# Area of Intersection of Polygons
#
# Algorithm based on http:#cap-lore.com/MathPhys/IP/
#
# Adapted 9-May-2006 by Lagado

import numpy


class Point:
    x = 0.0 
    y = 0.0
    def __init__(self,x,y):
        self.x = x
        self.y = y

class PolygonIntersect:
    class Box:
        min = None
        max = None
        def __init__(self,min,max):
            self.min = min
            self.max = max
            

    class Rng:
        mn = 0
        mx = 0
        def __init__(self,mn,mx):
            self.mn = mn
            self.mx = mx
            
    class IPoint:
        x = 0
        y = 0
            
    class Vertex:
        ip = None
        rx = None
        ry = None
        in_ = None
        

    gamut = 500000000.;
    mid = gamut / 2.;

    #--------------------------------------------------------------------------

    def range(self,points, bbox):
        for p in points:
            bbox.min.x = min(bbox.min.x, p.x)
            bbox.min.y = min(bbox.min.y, p.y)
            bbox.max.x = max(bbox.max.x, p.x)
            bbox.max.y = max(bbox.max.y, p.y)
            
        
    def area(self,a,p,q):
        A = p.x * q.y - p.y * q.x + a.x * (p.y - q.y) + a.y * (q.x - p.x)
        return A
    
    def ovl(self,p, q):
        return p.mn < q.mx and q.mn < p.mx

    #--------------------------------------------------------------------------

    ssss = 0;
    sclx = 0.0;
    scly = 0.0;

    def cntrib(self,f_x, f_y, t_x, t_y, w):
        self.ssss += w * (t_x - f_x) * (t_y + f_y) / 2

    def fit(self, x, fudge, B):
            cx = len(x)
            ix = list( range(cx+1) )
            c = cx;
            while (c > 0):
                c -= 1
                ix[c] = self.Vertex();
                ix[c].ip = self.IPoint();
                ix[c].ip.x = (int((x[c].x - B.min.x) * self.sclx - self.mid) & ~7) | fudge | (c & 1);
                ix[c].ip.y = (int((x[c].y - B.min.y) * self.scly - self.mid) & ~7) | fudge;

            ix[0].ip.y += cx & 1;
            ix[cx] = ix[0];

            c = cx;
            while (c > 0):
                c -= 1
                if ix[c].ip.x < ix[c + 1].ip.x:
                    ix[c].rx = self.Rng(ix[c].ip.x, ix[c + 1].ip.x)
                else: 
                    ix[c].rx = self.Rng(ix[c + 1].ip.x, ix[c].ip.x)

                if ix[c].ip.y < ix[c + 1].ip.y:
                    ix[c].ry = self.Rng(ix[c].ip.y, ix[c + 1].ip.y)
                else:
                    ix[c].ry = self.Rng(ix[c + 1].ip.y, ix[c].ip.y)
                ix[c].in_ = 0;
            return ix

    def cross(self,a, b, c, d, a1, a2, a3, a4):
            r1 = float(a1) / float( a1 + a2 )
            r2 = float(a3) / float( a3 + a4 )
            self.cntrib( (a.ip.x + r1 * (b.ip.x - a.ip.x)),\
                    (a.ip.y + r1 * (b.ip.y - a.ip.y)),\
                    b.ip.x, b.ip.y, 1)
            self.cntrib( d.ip.x, d.ip.y,\
                    (c.ip.x + r2 * (d.ip.x - c.ip.x)),\
                    (c.ip.y + r2 * (d.ip.y - c.ip.y)),\
                    1)
            a.in_ += 1
            c.in_ -= 1

    def inness(self,P, cP, Q, cQ):
        s = 0;
        c = cQ;
        p = P[0].ip;

        while ( c > 0):
            c -= 1
            if (Q[c].rx.mn < p.x and p.x < Q[c].rx.mx):
                sgn = 0 < self.area(p, Q[c].ip, Q[c + 1].ip)
                if (sgn != ( Q[c].ip.x < Q[c + 1].ip.x)):
                    s += 0
                elif sgn:
                    s += -1
                else:
                    s += 1

        for j in range (cP):
            if (s != 0):
                self.cntrib(P[j].ip.x, P[j].ip.y,
                        P[j + 1].ip.x, P[j + 1].ip.y, s);
            s += P[j].in_;

    #-------------------------------------------------------------------------

    def inter(self, a, b):
            na = len(a)
            nb = len(b)
            ipa = []
            ipb = []
            bbox = self.Box( Point(1e6, 1e6), Point(-1e6, -1e6) );

            if (na < 3 or nb < 3): return 0;

            self.range(a, bbox);
            self.range(b, bbox);

            rngx = bbox.max.x - bbox.min.x;
            self.sclx = self.gamut / rngx;
            rngy = bbox.max.y - bbox.min.y;
            self.scly = self.gamut / rngy;
            ascale = self.sclx * self.scly;

            ipa = self.fit(a, 0, bbox);
            ipb = self.fit(b, 2, bbox);

            for j in range(0,na):
                for k in range(0,nb):
                    if (self.ovl(ipa[j].rx, ipb[k].rx) and self.ovl(ipa[j].ry, ipb[k].ry)):
                        a1 = -self.area(ipa[j].ip, ipb[k].ip, ipb[k + 1].ip);
                        a2 =  self.area(ipa[j + 1].ip, ipb[k].ip, ipb[k + 1].ip);
                        o = a1 < 0;
                        if (o == (a2 < 0)):
                            a3 =  self.area(ipb[k].ip, ipa[j].ip, ipa[j + 1].ip);
                            a4 = -self.area(ipb[k + 1].ip, ipa[j].ip,ipa[j + 1].ip);
                            if ((a3 < 0) == (a4 < 0)):
                                if (o):
                                    self.cross(ipa[j], ipa[j + 1], ipb[k], ipb[k + 1],a1, a2, a3, a4);
                                else:
                                    self.cross(ipb[k], ipb[k + 1], ipa[j], ipa[j + 1],a3, a4, a1, a2);

            self.inness(ipa, na, ipb, nb);
            self.inness(ipb, nb, ipa, na);

            return self.ssss / ascale;

#  return the area of intersection of two polygons
#  Note: the area result has little more accuracy than a float
#  This is true even if the polygon is specified with doubles.
def intersectionArea(a, b):
    polygonIntersect = PolygonIntersect();
    return polygonIntersect.inter(a, b);

#--------------------------------------------------------------------------
#-------------------------------------------------------------------------
# test the code

def toPointsArray(inarray):
    A = []
    for a in inarray:
        A.append( Point( a[0], a[1] ) )
    return A;

def main():
    def trial(a, b):
        A = toPointsArray(a);
        B = toPointsArray(b);
        ab = intersectionArea(A, B)
        aa = intersectionArea(A, A)
        # print A[0].x, A[0].y
        print( "Result: ", ab, " ", aa )

    a1 = [[2,3], [2,3], [2,3], [2,4], [3,3], [2,3], [2,3]]
    b1 = [[1,1], [1,4], [4,4], [4,1], [1,1]] # 1/2, 1/2
    # The redundant vertices above are to provoke errors
    # as good test cases should.
    # It is not necessary to duplicate the first vertex at the end.
    a2 = [[1,7], [4,7], [4, 6], [2,6], [2, 3], [4,3], [4,2], [1,2]]
    b2 = [[3,1], [5,1], [5,4], [3,4], [3,5], [6,5], [6,0], [3,0]] # 0, 9
    a3 = [[1,1], [1,2], [2,1], [2,2]]
    b3 = [[0,0], [0,4], [4,4], [4,0]] # 0, 1/2
    a4 = [[0,0], [3,0], [3,2], [1,2], [1,1], [2,1], [2,3], [0,3]]
    b4 = [[0,0], [0,4], [4,4], [4,0]] # -9, 11
    a5 = [[0,0], [1,0], [0,1]]
    b5 = [[0,0], [0,1], [1,1], [1,0]] # -1/2, 1/2
    a6 = [[1, 3] , [2, 3] , [2, 0] , [1, 0] ]
    b6 = [[0, 1] , [3, 1] , [3, 2] , [0, 2] ] # -1, 3
    a7 = [[0,0], [0,2], [2,2], [2,0]]
    b7 = [[1, 1], [3, 1], [3, 3], [1, 3]] # -1, 4
    a8 = [[0,0], [0,4], [4,4], [4,0]]
    b8 = [[1,1], [1,2], [2,2], [2,1]] # 1, 16
    def plotpp(a,b):
        import pylab
        for i in range(len(a)-1):
            pylab.plot( [a[i][0],a[i+1][0]], [a[i][1],a[i+1][1]] ,'k-o')
        pylab.plot( [ a[i+1][0],a[0][0]],[a[i+1][1],a[0][1]] ,'k-o')
        for i in range(len(b)-1):
            pylab.plot( [b[i][0],b[i+1][0]], [b[i][1],b[i+1][1]] ,'b-o')
        pylab.plot( [ b[i+1][0],b[0][0]],[b[i+1][1],b[0][1]] ,'b-o')
        pylab.show()
    trial(a1, b1)
    trial(a2, b2)
    trial(a3, b3)
    trial(a4, b4)
    trial(a5, b5)
    trial(a6, b6)
    trial(a7, b7)
    trial(a8, b8)
    #plotpp(a5,b5)

if __name__ == '__main__': 
        main()
