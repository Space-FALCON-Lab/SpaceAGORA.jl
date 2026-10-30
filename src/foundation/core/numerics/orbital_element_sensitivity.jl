module OrbitalElementSensitivity
using LinearAlgebra
export orbital_element_jacobian

skew(x)=[0. -x[3] x[2];x[3] 0. -x[1];-x[2] x[1] 0.]

"""Analytical d[a,e,i,RAAN,argument of periapsis]/d[r,v], SI and radians.

Adapted from the paper ZIP's classical_element_sensitivity.jl. The classical
angular derivatives require a noncircular, nonequatorial orbit.
"""
function orbital_element_jacobian(r,v,mu)
    radius=norm(r);h=cross(r,v);hm=norm(h);node=cross([0.,0.,1.],h);nm=norm(node)
    ev=cross(v,h)/mu-r/radius;e=norm(ev)
    e>1e-10 && nm/hm>1e-10 || throw(ArgumentError("Classical element derivatives require noncircular, nonequatorial orbits."))
    I3=Matrix{Float64}(I,3,3);dr=hcat(I3,zeros(3,3));dv=hcat(zeros(3,3),I3)
    dh=-skew(v)*dr+skew(r)*dv;dn=skew([0.,0.,1.])*dh
    de=(-skew(h)*dv+skew(v)*dh)/mu-(I3/radius-r*r'/radius^3)*dr
    hh=h/hm;nn=node/nm;ee=ev/e
    dhh=(I3-hh*hh')*dh/hm;dnn=(I3-nn*nn')*dn/nm;dee=(I3-ee*ee')*de/e
    a=inv(2/radius-dot(v,v)/mu)
    da=2a^2/mu*(mu*r'/radius^3*dr+v'*dv)
    di=(h[3]*(h[1]*dh[1:1,:]+h[2]*dh[2:2,:])/nm-nm*dh[3:3,:])/hm^2
    dnode=(node[1]*dn[2:2,:]-node[2]*dn[1:1,:])/nm^2
    x=dot(nn,ee);y=dot(hh,cross(nn,ee))
    dx=ee'*dnn+nn'*dee
    dy=cross(nn,ee)'*dhh+hh'*(-skew(ee)*dnn+skew(nn)*dee)
    return vcat(da,ee'*de,di,dnode,(x*dy-y*dx)/(x^2+y^2))
end
end
