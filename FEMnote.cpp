vega_FEM_initiate_cpp(void)
    initConfigurations()
        // read vega.config
        // volumetircMeshFilename = mesh.veg
        // solver = implicitNewmark
        // bodyType = 2
    initSimulation()
        // deformableObject = NULL, 0: STVK, 1: COROTLINFEM, 2: LINFEM, 3: MESSSRPING 4: (INVERTIBLEFEM), 5: UNSPECIFIED
        // invertibleMaterial = 0: (INV_STVK), 1: INV_NEOHOOKEAN, 2: MooneyRivlin
        // massSprintSystemSource = (NULL), 0: OBJ, 1: TETMESH, 2: CUBICMESH, 3: CHAIN
        // read mesh.veg
            // volumetricMesh: 
            // numVertices = 348, numElements = 801
            // vertices[0:347]: (x,y,z), elements[0:800]: (v1,v2,v3,v4)
        // n = 348
        // set meshGraph, LaplacianDampingMatrix, dampingLaplacianCoef
        computeMassMatrix(volumetricMesh, &massMatrix, inflate3Dim)
        stVKInternalForces = new StVKInternalForces(volumetricMesh, precomputedIntegrals, addGravity, g)
        stVKStiffnessMatrix = new StVKStiffnessMatrix(stVKInternalForces)
        fixedVertices[0:nFixed] // read mesh.bou
        fixedDOFS[0:3nFixed-1]  // index # of fixed dof
        allocate u, u_base, uvel, uvel_base, uaccel, uaccel_base
        allocate f_ext, f_extBase, f_int, f_int_base, verticesRelation_vega_IB
        tetMesh = dynamic_cast<TetMesh*>(volumetricMesh);
        isotropicMaterial = new StVKIsotropicMaterial(tetMesh, enableCompressionResistance, compressionResistance)
        isotropicHyperelasticFEM = new IsotropicHyperelasticFEM(tetMesh, isotropicMaterial, inversionThreshold, addGravity, g)
        implicitNewmarkSparse = new ImplicitNewmarkSparse(3*n, timeStep, massMatrix, forceModel, positiveDefinite, numFixedDOFs, fixedDOFs,
            dampingMassCoef, dampingStiffnessCoef, maxIterations, epsilon, newmarkBeta, newmarkGamma, numSolverThreads);
        integratorBase = implicitNewmarkSparse
        integratorBaseSparse->SetDampingMatrix(LaplacianDampingMatrix);
        integratorBase->ResetToRest();
        integratorBase->SetState(uInitial, velInitial);
        integratorBase->SetTimestep(timeStep/substepsPerTimeStep);
    // write TET_ORIGINAL.obj
    readFlowGridXYZ()
        // read xgrid.dat, ygrid.dat to flowGrid[0:2][0:nxc-1] (Float)
    readSurfaceTriangleMesh()
        // read vegaSurfaceMesh.dat to surfaceTriangleMesh[0:3nsTri-1] (Int)
    // nVertices = 174
vega_interpolateIndexRatio_cpp(int *markerInterpolateIndex_,double *markerInterpolateRatio_)
    obtainBodyNormalVectorCW()
        // bodyNormalVectorCW[0:3nsVer-1] (Float), norm vector of the surface
    obtainMarkerInterpolateIndexRatio(markerInterpolateIndex_,markerInterpolateRatio_)
vega_deformation_cpp(double *markerPressure_,double *markerInterpolateVelocity_,double *bodyMotion_)
    obtainMarkerForceFromPressureAndVelocity(markerPressure_,markerInterpolateVelocity_,bodyMotion_);
        // markerForces[0:3nsVer-1] (Float), pressure + shear
        // save markerForce to Vega_Force_total.dat
    // f_ext[0:3n]
    // integratorBase
        // heavingAmplitude = 0.05, heavingFrequency = 1.5
        // externalForces
        SetqState()
            // q, qvel, qaccel = u_base, uvel_base, uaacel_base
        DoTimestep()
        GetqState()
            // u, uvel, uaccel = q, qvel, qaccel
    // bodyMotion_[0:nsVel-1] (Float), surface velocity
vega_reNewBodyPosition_cpp()
    // integratorBase
        GetqState()
            // u_base, uvel_base, uaccel_base = q, qvel, qaccel
    // timestepCounter++




        



DoTimeStep()
    newmarkBeta = 0.25
    newmarkGamma = 0.5
    alpha1 = 1.0 / (NewmarkBeta * timestep * timestep)
    alpha2 = 1.0 / (NewmarkBeta * timestep)
    alpha3 = (1.0 - 2.0 * NewmarkBeta) / (2.0 * NewmarkBeta)
    alpha4 = NewmarkGamma / (NewmarkBeta * timestep)
    alpha5 = 1 - NewmarkGamma/NewmarkBeta
    alpha6 = (1.0 - NewmarkGamma / (2.0 * NewmarkBeta)) * timestep

    // set the initial guess of qaccel, qvel
    for(int i=0; i<3n; i++)
    {
        q_1[i] = q[i];
        qvel_1[i] = qvel[i];
        qaccel_1[i] = qaccel[i];
        qaccel[i] = alpha1 * (q[i] - q_1[i]) - alpha2 * qvel_1[i] - alpha3 * qaccel_1[i];
        qvel[i] = alpha4 * (q[i] - q_1[i]) + alpha5 * qvel_1[i] + alpha6 * qaccel_1[i];
    }
