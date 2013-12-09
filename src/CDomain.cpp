#include "CDomain.h"
#include <math.h>     //ds fabs, etc.
#include "writepng.h" //ds p2m png
#include <mpi.h>      //ds MPI
#include <iostream>



#define MPI_WORKTAG 1
#define MPI_DIETAG 2

//ds speed
static const double M_PI_SQUARED = M_PI*M_PI;

namespace Diffusion
{

//ds ctor/dtor
CDomain::CDomain( const double& p_dDiffusionCoefficient,
             const std::pair< double, double >& p_prBoundaries,
             const double& p_dBoundarySize,
             const unsigned int& p_uNumberOfGridPoints1D,
             const unsigned int& p_uNumberOfParticles,
             const double& p_dGridPointSpacing,
             const double& p_dTimeStepSize,
             const unsigned int& p_uRank,
             const unsigned int& p_uNumberOfTasks ) : m_dDiffusionCoefficient( p_dDiffusionCoefficient ),
                                                      m_prBoundaries( p_prBoundaries ),
                                                      m_dBoundarySize( p_dBoundarySize ),
                                                      m_uNumberOfGridPoints1D( p_uNumberOfGridPoints1D ),
                                                      m_uNumberOfParticles( p_uNumberOfParticles ),
                                                      m_dGridPointSpacing( p_dGridPointSpacing ),
                                                      m_dTimeStepSize( p_dTimeStepSize ),
                                                      m_dEpsilon( 2*p_dGridPointSpacing ),
                                                      m_dVolP( p_dGridPointSpacing*p_dGridPointSpacing ),
                                                      m_PSEFactor( p_dTimeStepSize*p_dDiffusionCoefficient/( m_dEpsilon*m_dEpsilon )*m_dVolP ),
                                                      m_strLogHeatDistribution( "" ),
                                                      m_strLogNorms( "" ),
                                                      m_uNumberOfTasks( p_uNumberOfTasks ),
                                                      m_uRank( p_uRank ),
                                                      m_uLoopSize( p_uNumberOfGridPoints1D/( p_uNumberOfTasks-1 ) )

{
    //ds only allocate grid if we are the master
    if( 0 == m_uRank )
    {
        //ds allocate memory for the data structure
        m_gridHeat = new double*[m_uNumberOfGridPoints1D];

        //ds for each element
        for( unsigned int u = 0; u < m_uNumberOfGridPoints1D; ++u )
        {
            m_gridHeat[u] = new double[m_uNumberOfGridPoints1D];
        }

        //ds initialize grid
        setInitialHeatDistribution( );
    }
};

CDomain::~CDomain( )
{
    //ds only for the master
    if( 0 == m_uRank )
    {
        //ds deallocate heat structure
        for( unsigned int u = 0; u < m_uNumberOfGridPoints1D; ++u )
        {
            delete[] m_gridHeat[u];
        }

        delete[] m_gridHeat;
    }
};

//ds accessors
void CDomain::updateHeatDistributionNumerical( )
{
    //ds heat change for current time step
    double gridHeatChangePSE[m_uNumberOfGridPoints1D][m_uNumberOfGridPoints1D];

    //ds for all grid points
    for( unsigned int uk = 0; uk < m_uNumberOfGridPoints1D; ++uk )
    {
        for( unsigned int vk = 0; vk < m_uNumberOfGridPoints1D; ++vk )
        {
            //ds get 2d vector of current coordinates
            const double dXk[2] = { uk*m_dGridPointSpacing, vk*m_dGridPointSpacing };

            //ds inner sum of formula
            double dInnerSum( 0.0 );

            //ds loop 20x20
            for( int i = -10; i <= 10; ++i )
            {
                for( int j = -10; j <= 10; ++j )
                {
                    //ds get current indexes up, vp
                    int up = uk + i;
                    int vp = vk + j;

                    //ds if we are not ourself (unsigned overflow no problem here)
                    if( uk != static_cast< unsigned int >( up ) && vk != static_cast< unsigned int >( vp ) )
                    {
                        //ds offset value (required for spaced positions)
                        double dOffsetU = 0.0;
                        double dOffsetV = 0.0;

                        //ds check boundary
                        if( 0 > up                        ){ up += m_uNumberOfGridPoints1D; dOffsetU = -m_dBoundarySize; } //ds moving to negative boundary
                   else if( m_uNumberOfGridPoints1D <= up ){ up -= m_uNumberOfGridPoints1D; dOffsetU = m_dBoundarySize;  } //ds moving to positive boundary
                        if( 0 > vp                        ){ vp += m_uNumberOfGridPoints1D; dOffsetV = -m_dBoundarySize; } //ds moving to negative boundary
                   else if( m_uNumberOfGridPoints1D <= vp ){ vp -= m_uNumberOfGridPoints1D; dOffsetV = m_dBoundarySize;  } //ds moving to positive boundary

                        //ds get 2d vector of current coordinates
                        const double dXp[2] = { ( up*m_dGridPointSpacing + dOffsetU ), ( vp*m_dGridPointSpacing + dOffsetV ) };

                        //ds compute inner sum
                        dInnerSum += ( m_gridHeat[up][vp] - m_gridHeat[uk][vk] )*getKernelEta( dXp, dXk );
                    }
                }
            }

            //ds add final part of formula and save in temporary grid
            gridHeatChangePSE[uk][vk] = m_PSEFactor*dInnerSum;
        }
    }

    //ds copy all computed values to original grid
    for( unsigned int u = 0; u < m_uNumberOfGridPoints1D; ++u )
    {
        for( unsigned int v = 0; v < m_uNumberOfGridPoints1D; ++v )
        {
            m_gridHeat[u][v] += gridHeatChangePSE[u][v];
        }
    }
}

void CDomain::updateHeatDistributionNumericalMASTER( )
{
    //ds get portions of the data to send to the slaves
    std::vector< double** > vecGridUnits;

    //ds fill the vector
    for( unsigned int u = 1; u < m_uNumberOfTasks; ++u )
    {
        //ds get a copy of our heat grid - TODO: reduce massive memory expenses
        vecGridUnits.push_back( getCopyOfHeatGrid( ) );
    }

    //ds send the data to the workers
    for( int iRank = 1; iRank < m_uNumberOfTasks; ++iRank )
    {
        int* currentGrid = new int[10];

        std::cout << "current grid: " << currentGrid[0] << std::endl;
        std::cout << "current grid: " << currentGrid[10-1] << std::endl;

        std::cout << "sending to: " << iRank << std::endl;

        //ds send respective grid unit to slave
        MPI_Send( &currentGrid, 10, MPI_INT, iRank, MPI_WORKTAG, MPI_COMM_WORLD );
    }

    while( true ){ };

    //ds wait for all results from workers
    for( int iRank = 1; iRank < m_uNumberOfTasks; ++iRank )
    {
        //ds result grid from worker
        double** gridHeadResult( 0 );

        //ds get the result
        MPI_Recv( &(gridHeadResult[0][0]), m_uNumberOfGridPoints1D*m_uNumberOfGridPoints1D, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, 0 );

        //ds update the main grid
        for( unsigned int u = 0; u < m_uNumberOfGridPoints1D; ++u )
        {
            for( unsigned int v = 0; v < m_uNumberOfGridPoints1D; ++v )
            {
                //ds this is possible because slaves will have 0 values for parts they do not contribute to
                m_gridHeat[u][v] += gridHeadResult[u][v];
            }
        }

        //ds release inner elements
        for( unsigned int u = 0; u < m_uNumberOfGridPoints1D; ++u )
        {
            delete[] gridHeadResult[u];
        }

        //ds and main element
        delete[] gridHeadResult;
    }

    //ds tell the slaves to exit
    for( int iRank = 1; iRank < m_uNumberOfTasks; ++iRank )
    {
        //ds die command
        MPI_Send( 0, 0, MPI_INT, iRank, MPI_DIETAG, MPI_COMM_WORLD );
    }

    //ds done
    std::cout << "master complete" << std::endl;
    getchar( );
}

void CDomain::updateHeatDistributionNumericalSLAVE( )
{
    //ds status
    MPI_Status mpiStatus;

    //ds get looping range
    const unsigned int uIndexStart( ( m_uRank-1 )*m_uLoopSize );
    unsigned int uIndexEnd( m_uRank*m_uLoopSize - 1 );

    //ds if we are the final slave the index end has to be max (sometimes one can enter odd numbers of processes/particles so they can not be divided properly)
    if( m_uNumberOfTasks == m_uRank + 1 )
    {
        //ds set it to the end
        uIndexEnd = m_uNumberOfGridPoints1D;
    }

    //ds until master stops us
    while( true )
    {
        //ds heat grid to work with
        int* gridHeat;

        std::cout << "receiving: " << m_uRank << std::endl;

        //ds receive message from the master
        MPI_Recv( &gridHeat, 10, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &mpiStatus );

        std::cout << "task: " << m_uRank << " received heat grid: " << uIndexStart << " to " << uIndexEnd << std::endl;

        std::cout << "gid[0][0]: " << std::endl;
        std::cout << gridHeat[0] << std::endl;

        std::cout << "over" << std::endl;

        //ds if its the stop tag terminate
        if( MPI_DIETAG == mpiStatus.MPI_TAG )
        {
            //ds escape loop
            return;
        }

        //ds heat change for current time step
        double gridHeatChangePSE[m_uNumberOfGridPoints1D][m_uNumberOfGridPoints1D];

        //ds for all grid points in our range
        for( unsigned int uk = uIndexStart; uk < uIndexEnd; ++uk )
        {
            for( unsigned int vk = uIndexStart; vk < uIndexEnd; ++vk )
            {
                //ds get 2d vector of current coordinates
                const double dXk[2] = { uk*m_dGridPointSpacing, vk*m_dGridPointSpacing };

                //ds inner sum of formula
                double dInnerSum( 0.0 );

                //ds loop 20x20
                for( int i = -10; i <= 10; ++i )
                {
                    for( int j = -10; j <= 10; ++j )
                    {
                        //ds get current indexes up, vp
                        int up = uk + i;
                        int vp = vk + j;

                        //ds if we are not ourself (unsigned overflow no problem here)
                        if( uk != static_cast< unsigned int >( up ) && vk != static_cast< unsigned int >( vp ) )
                        {
                            //ds offset value (required for spaced positions)
                            double dOffsetU = 0.0;
                            double dOffsetV = 0.0;

                            //ds check boundary
                            if( 0 > up                        ){ up += m_uNumberOfGridPoints1D; dOffsetU = -m_dBoundarySize; } //ds moving to negative boundary
                       else if( m_uNumberOfGridPoints1D <= up ){ up -= m_uNumberOfGridPoints1D; dOffsetU = m_dBoundarySize;  } //ds moving to positive boundary
                            if( 0 > vp                        ){ vp += m_uNumberOfGridPoints1D; dOffsetV = -m_dBoundarySize; } //ds moving to negative boundary
                       else if( m_uNumberOfGridPoints1D <= vp ){ vp -= m_uNumberOfGridPoints1D; dOffsetV = m_dBoundarySize;  } //ds moving to positive boundary

                            //ds get 2d vector of current coordinates
                            const double dXp[2] = { ( up*m_dGridPointSpacing + dOffsetU ), ( vp*m_dGridPointSpacing + dOffsetV ) };

                            //ds compute inner sum
                            //dInnerSum += ( gridHeat[up][vp] - gridHeat[uk][vk] )*getKernelEta( dXp, dXk );
                        }
                    }
                }

                //ds add final part of formula and save in temporary grid
                //gridHeatChangePSE[uk][vk] = m_PSEFactor*dInnerSum + gridHeat[uk][vk];
            }
        }

        //ds copy all computed values to original grid
        for( unsigned int u = uIndexStart; u < uIndexEnd; ++u )
        {
            for( unsigned int v = uIndexStart; v < uIndexEnd; ++v )
            {
                //ds this assignment causes zero values for no changes
                //gridHeat[u][v] = gridHeatChangePSE[u][v];
            }
        }

        //ds send the result back
        MPI_Send( &gridHeat, m_uNumberOfGridPoints1D*m_uNumberOfGridPoints1D, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    }
}

void CDomain::updateHeatDistributionAnalytical( const double& p_dCurrentTime )
{
    //ds for all grid points
    for( unsigned int u = 0; u < m_uNumberOfGridPoints1D; ++u )
    {
        for( unsigned int v = 0; v < m_uNumberOfGridPoints1D; ++v )
        {
            //ds get the exact heat at this point
            m_gridHeat[u][v] = getHeatAnalytical( u*m_dGridPointSpacing, v*m_dGridPointSpacing, p_dCurrentTime );
        }
    }
}

void CDomain::saveHeatGridToStream( )
{
    //ds buffer for snprintf
    char chBuffer[16];

    //ds add each element
    for( unsigned int u = 0; u < m_uNumberOfGridPoints1D; ++u )
    {
        for( unsigned int v = 0; v < m_uNumberOfGridPoints1D; ++v )
        {
            //ds get the integrals stream
            std::snprintf( chBuffer, 16, "%f", m_gridHeat[u][v] );

            //ds add buffer and space to string
            m_strLogHeatDistribution += chBuffer;
            m_strLogHeatDistribution += " ";
        }

        //ds add new line for new row
        m_strLogHeatDistribution += '\n';
    }
}

void CDomain::saveNormsToStream( const double& p_dCurrentTime )
{
    //ds get total heat
    double dTotalHeat( 0.0 );

    //ds norms
    double dLInfinity( 0.0 );
    double dL1( 0.0 );
    double dL2( 0.0 );

    //ds for all grid points
    for( unsigned int u = 0; u < m_uNumberOfGridPoints1D; ++u )
    {
        for( unsigned int v = 0; v < m_uNumberOfGridPoints1D; ++v )
        {
            //ds add the current heat value
            dTotalHeat += m_gridHeat[u][v];

            //ds calculate the error (numerical - analytical) - the absolute value is already here taken since it does not change L
            const double dErrorAbsolute( fabs( m_gridHeat[u][v] - getHeatAnalytical( u*m_dGridPointSpacing, v*m_dGridPointSpacing, p_dCurrentTime ) ) );

            //ds check for LInf
            if( dLInfinity < dErrorAbsolute )
            {
                dLInfinity = dErrorAbsolute;
            }

            //ds L1, L2
            dL1 += dErrorAbsolute;
            dL2 += dErrorAbsolute*dErrorAbsolute;
        }
    }

    //ds scale L1, L2
    dL1 /= m_uNumberOfParticles;
    dL2 /= m_uNumberOfParticles;
    dL2 = sqrt( dL2 );

    //ds buffer for snprintf
    char chBuffer[64];

    //ds get the norms stream: E Linf L1 L2
    std::snprintf( chBuffer, 64, "%f %f %f %f", dTotalHeat, dLInfinity, dL1, dL2 );

    //ds append the buffer to our string
    m_strLogNorms += chBuffer;
    m_strLogNorms += "\n";
}

void CDomain::saveMeshToPNG( const unsigned int& p_dCurrentTimeStep, const unsigned int& p_uRate )
{
    //ds if the current step matches the rate
    if( 0 == p_dCurrentTimeStep%p_uRate )
    {
        //ds get mesh size
        const unsigned int uMeshSize( 2*m_uNumberOfGridPoints1D );

        //ds get the data
        unsigned char* chMesh( getP2M( uMeshSize ) );

        //ds construct picture name - buffer for snprintf
        char chBuffer[64];

        //ds format: image_number.png
        std::snprintf( chBuffer, 64, "bin/image_%.4u.png", 0 );

        //ds create png
        writePNG( chBuffer, uMeshSize, uMeshSize, chMesh );

        //ds free data
        delete chMesh;
    }
}

void CDomain::writeHeatGridToFile( const std::string& p_strFilename, const unsigned int& p_uNumberOfTimeSteps ) const
{
    //ds ofstream object
    std::ofstream ofsFile;

    //ds open the file for writing
    ofsFile.open( p_strFilename.c_str( ), std::ofstream::out );

    //ds if it worked
    if( ofsFile.is_open( ) )
    {
        //ds first dump setup information number of points and time steps
        ofsFile << m_uNumberOfGridPoints1D << " " << p_uNumberOfTimeSteps << "\n" << m_strLogHeatDistribution;
    }

    //ds close the file
    ofsFile.close( );
}

void CDomain::writeNormsToFile( const std::string& p_strFilename, const unsigned int& p_uNumberOfTimeSteps, const double& p_dTimeStepSize ) const
{
    //ds ofstream object
    std::ofstream ofsFile;

    //ds open the file for writing
    ofsFile.open( p_strFilename.c_str( ), std::ofstream::out );

    //ds if it worked
    if( ofsFile.is_open( ) )
    {
        //ds dump information to file
        ofsFile << p_uNumberOfTimeSteps << " " << p_dTimeStepSize << "\n" << m_strLogNorms;
    }

    //ds close the file
    ofsFile.close( );
}

//ds helpers
void CDomain::setInitialHeatDistribution( )
{
    //ds loop over all indexi
    for( unsigned int u = 0; u < m_uNumberOfGridPoints1D; ++u )
    {
        for( unsigned int v = 0; v < m_uNumberOfGridPoints1D; ++v )
        {
            //ds set initial heat value
            m_gridHeat[u][v] = sin( u*m_dGridPointSpacing*2*M_PI )*sin( v*m_dGridPointSpacing*2*M_PI );
        }
    }
}

double CDomain::getHeatAnalytical( const double& p_dX, const double& p_dY, const double& p_dT ) const
{
    //ds formula
    return sin( p_dX*2*M_PI )*sin( p_dY*2*M_PI )*exp( -8*m_dDiffusionCoefficient*M_PI_SQUARED*p_dT );
}

double CDomain::getKernelEta( const double p_dXp[2], const double p_dXk[2] ) const
{
    //ds compute distance (using pow for readability)
    const double dDistance = sqrt( ( p_dXp[0] - p_dXk[0] )*( p_dXp[0] - p_dXk[0] ) + ( p_dXp[1] - p_dXk[1] )*( p_dXp[1] - p_dXk[1] ) );

    //ds if we are out of the cutoff
    if( 5*m_dEpsilon < dDistance )
    {
        return 0.0;
    }
    else
    {
        //return kernel function
        return 16.0/( M_PI_SQUARED*( pow( dDistance, 8 ) + 1.0 ) );
    }
}

double CDomain::getKernelW( const double& p_dLambda ) const
{
    //ds get absolute lambda
    const double dLambdaAbs = fabs( p_dLambda );

    //ds cases
    if( 0 <= dLambdaAbs && dLambdaAbs < 1 )
    {
        return ( 1.0 - 5.0/2.0*dLambdaAbs*dLambdaAbs + 3.0/2.0*dLambdaAbs*dLambdaAbs*dLambdaAbs );
    }
    else if( 1 <= dLambdaAbs && dLambdaAbs < 2 )
    {
        return ( 1.0/2.0*( 2.0 - dLambdaAbs )*( 2.0 - dLambdaAbs )*( 1.0 - dLambdaAbs ) );
    }
    else
    {
        return 0.0;
    }
}

double CDomain::clamp(const double &x, const double &a, const double &b) const
{
    return x < a ? a : x > b ? b : x;
}

unsigned char* CDomain::getP2M( const unsigned int& p_uMeshSize ) const
{
    //ds h value
    const double dH( 1.0/p_uMeshSize );

    //ds allocate data array
    unsigned char* chMesh = new unsigned char[p_uMeshSize*p_uMeshSize*4];

    //ds set data - for each data element
    for( unsigned int i = 0; i < p_uMeshSize; ++i )
    {
        for( unsigned int j = 0; j < p_uMeshSize; ++j )
        {
            //ds index in mesh
            unsigned int uIndexMesh = 4.0*( i*p_uMeshSize + j );

            //ds concentration on mesh
            double dConcentration( 0.0 );

            //ds for all particles
            for( unsigned int u = 0; u < m_uNumberOfGridPoints1D; ++u )
            {
                for( unsigned int v = 0; v < m_uNumberOfGridPoints1D; ++v )
                {
                    //ds get lambdas
                    const double dLambdai = u*m_dGridPointSpacing/dH - i;
                    const double dLambdaj = v*m_dGridPointSpacing/dH - j;

                    //ds compute concentration
                    dConcentration += m_gridHeat[u][v]*getKernelW( dLambdai )*getKernelW( dLambdaj );
                }
            }

            //ds set the mesh values - SNIPPET
            chMesh[uIndexMesh+0] = clamp( gamma( dConcentration ), 0.0, 0.1 )*0xFFu;  // red
            chMesh[uIndexMesh+1] = clamp( 1.0, 0.0, 1.0 )*0xFFu;                      // green
            chMesh[uIndexMesh+2] = clamp( gamma( -dConcentration ), 0.0, 0.1 )*0xFFu; // blue
            chMesh[uIndexMesh+3] = 0xFFu;                                             // alpha (opacity, 0xFFu = 255)
        }
    }

    return chMesh;
}

double** CDomain::getCopyOfHeatGrid( ) const
{
    //ds allocate memory for the data structure
    double** gridHeat = new double*[m_uNumberOfGridPoints1D];

    //ds for each element
    for( unsigned int u = 0; u < m_uNumberOfGridPoints1D; ++u )
    {
        gridHeat[u] = new double[m_uNumberOfGridPoints1D];
    }

    //ds take over all values
    for( unsigned int u = 0; u < m_uNumberOfGridPoints1D; ++u )
    {
        for( unsigned int v = 0; v < m_uNumberOfGridPoints1D; ++v )
        {
            //ds get value from real grid
            gridHeat[u][v] = m_gridHeat[u][v];
        }
    }

    //ds return the copy
    return gridHeat;
}

} //namespace Diffusion
