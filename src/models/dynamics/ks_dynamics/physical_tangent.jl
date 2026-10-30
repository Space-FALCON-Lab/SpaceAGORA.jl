"""Fictitious-time KS input matrix for inertial acceleration (m/s^2)."""
function ks_acceleration_input_matrix(state)
    p=state[1:4];q=state[5:8];rho=dot(p,p)
    B=zeros(10,3)
    B[5:8,:].=0.5rho*transpose(ks_lambda_matrix(p))
    B[9,:].=-rho*ks_velocity(p,q)
    return B
end

"""Physical-time tangent and analytical constrained lift, with h=-specific energy.

The nominal force model is the existing point-mass plus J2 KS model. Input
columns are arbitrary inertial accelerations. This is not an active-link tangent.
"""
function ks_physical_tangent(state,params)
    evaluation=evaluate_ks_dynamics_and_jacobians(state,params;compute_input_jacobian=Val(false))
    fs=Vector(evaluation.derivative);As=Matrix(evaluation.state)
    p=Vector(state[1:4]);q=Vector(state[5:8]);rho=dot(p,p)
    r=ks_position(p);v=ks_velocity(p,q)
    pr,vp,vq=ks_kinematics_jacobians(p,q)[1:3]
    H=zeros(6,10);H[1:3,1:4].=pr;H[4:6,1:4].=vp;H[4:6,5:8].=vq
    ft=fs/rho;At=As/rho;At[:,1:4].-=fs*(2p)'/rho^2
    Bt=ks_acceleration_input_matrix(state)/rho
    pd=ft[1:4];qd=ft[5:8];rhod=2dot(p,pd);acc=(H*ft)[4:6]
    Hd=zeros(6,10)
    Hd[1:3,1:4].=2ks_lambda_matrix(pd)
    Hd[4:6,1:4].=2ks_lambda_matrix(qd)/rho-2rhod*ks_lambda_matrix(q)/rho^2-
        2(acc*p'+v*pd')/rho+2rhod*v*p'/rho^2
    Hd[4:6,5:8].=2ks_lambda_matrix(pd)/rho-2rhod*ks_lambda_matrix(p)/rho^2
    G=zeros(10,6);G[1:4,1:3].=ks_lambda_matrix(p)'/(2rho)
    for k in 1:3
        G[5:8,k].=0.5ks_lambda_matrix(G[1:4,k])'*v
    end
    G[5:8,4:6].=0.5ks_lambda_matrix(p)'
    G[9,:].=vcat(-params.μ*r/rho^3,-v)
    return (;A_time=At,B_time=Bt,H,Hdot=Hd,G,
        A_cartesian=(Hd+H*At)*G,B_cartesian=H*Bt,acceleration=acc)
end
