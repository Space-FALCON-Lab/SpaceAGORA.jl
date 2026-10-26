using LoopVectorization
using LinearAlgebra
using SatelliteToolboxTransformations: ned_to_ecef
using SatelliteToolboxGeomagneticField
using CSV
using DataFrames
using SpecialFunctions: loggamma
using ...QuaternionMath
# import .config
# Define delta function
δ(x,y) = ==(x,y)
δ(x) = δ(x,0)
