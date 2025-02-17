# ABTS.jl

# Setting up the Docker Environment (Currently not working for Apple Silicon Macs)
To ensure that all the package versions are consistent, a Docker environment has been configured for use with ABTS.jl. This guarantees consistent results between computers. This is needed for running on Apple Silicon Macs because of how GRAM was compiled. For other virtualization options, see [GRAM Setup](https://github.com/Space-FALCON-Lab/Aerobraking-Trajectory-Simulator/tree/GRAM-updates?tab=readme-ov-file#gram-setup). The process of setting it up is as follows:
1. Download and sign in to Docker Desktop ([link](https://www.docker.com/get-started/)). If on Linux, follow these instructions to sign in ([link](https://docs.docker.com/desktop/setup/sign-in/))
2. In VSCode, download the Docker and Dev Container extensions
3. Open the ABTS.jl directory in VSCode and click on the remote window icon in the bottom left corner
4. Select "Reopen in container". This will open a window with the ABTS.jl repository open in the preconfigured Docker environment.
5. Open a terminal window in VSCode. Start Julia, enter "]" to activate the package environment, and enter ```activate .ABTS```, then ```instantiate```. This will activate the Julia project environment that has been set up and install all the packages needed to run ABTS.jl
6. Now, you should be able to run any existing scenario or create and run a new scenario