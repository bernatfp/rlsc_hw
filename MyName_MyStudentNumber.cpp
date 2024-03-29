/*
 * BaxterClient.cpp
 *
 *  Created on: 24 Jan 2014
 *      Author: Vladimir Ivan
 */

// Make sure to have the server side running in V-REP!

#include "BaxterTools.h"

int main(int argc,char* argv[])
{
  std::cout << "Starting..." << std::endl;
  // Create the robot interface object
  BaxterTools bax;
  // Connect to the simulator
  if(argc==2)
  {
    bax.Connect(argv[1]);
  }
  else
  {
    bax.Connect("localhost");
  }

  std::cout << "Connected" << std::endl;

  // Start the simulation
	bax.StartSimulation();

  std::cout << "Simulation started" << std::endl;

	Eigen::VectorXd q = Eigen::VectorXd::Zero(18); // Joint angles
	Eigen::VectorXd x; // End-effector position
	Eigen::VectorXd target; // Target positions
	Eigen::MatrixXd J; // Jacobian matrix


  std::cout << "Vars declared" << std::endl;

	//////////////////////////////////////////////////////////////////////
	// Constants for homework
	Eigen::MatrixXd Winv = Eigen::MatrixXd::Zero(7,7); // Weighting matrix
	for(int i=0;i<7;i++) Winv(i,i) = ((double)i)/6.0+0.1;

	Eigen::MatrixXd C = Eigen::MatrixXd::Identity(3,3)*1e3; // Regularisation
    Eigen::MatrixXd Cinv = Eigen::MatrixXd::Identity(3,3)*1e-3;

	Eigen::VectorXd qstart1(18); // Starting pose 1
	Eigen::VectorXd qstart2(18); // Starting pose 2
	Eigen::VectorXd qstart3(18); // Starting pose 3
	qstart1 << M_PI/4.0,0,0,M_PI/2.0,0,0,0,0,0,      -M_PI/4.0,0,0,M_PI/2.0,0,0,0,0,0;
  qstart2 << -M_PI/4.0,0,0,M_PI/2.0,M_PI/2.0,M_PI/2.0,0,0,0,   M_PI/4.0,0,0,M_PI/2.0,-M_PI/2.0,M_PI/2.0,0,0,0;
  qstart3 << M_PI/4.0,-M_PI/2.0,0,M_PI/2.0,-M_PI/4.0,-M_PI/4.0,0,0,0,      -M_PI/4.0,-M_PI/2.0,0,M_PI/2.0,M_PI/4.0,-M_PI/4.0,0,0,0;

	Eigen::VectorXd q_comf1(18); // Comfortable pose 1
	Eigen::VectorXd q_comf2(18); // Comfortable pose 2
	q_comf1 << 0,-0.5,0,0.5,0,1.5,0,0,0,  0,-0.5,0,0.5,0,1.5,0,0,0 ;
	q_comf2 << -20.0/180.0*M_PI, 40.0/180.0*M_PI, 70.0/180.0*M_PI, 90.0/180.0*M_PI, 0,0.5,0,0,0, 20.0/180.0*M_PI, 40.0/180.0*M_PI, -70.0/180.0*M_PI, 90.0/180.0*M_PI, 0,0.5,0,0,0;

	//////////////////////////////////////////////////////////////////////


  std::cout << "Constants filled" << std::endl;

	// Loop until 'q' gets pressed
	char key=0;
	while(key!='q')
  {
    std::cout << "Press key" << std::endl;
  	// Get pressed key
  	key=bax.GetKey();
    std::cout << "Key extracted: " << key << std::endl;

    // Put your code here //////////////////////////

    //Part A

    std::cout << "Part A" << std::endl;

    //Init variables
    int idx; //aux index
    Eigen::MatrixXd costs = Eigen::MatrixXd::Zero(8,4); //Costs
    Eigen::VectorXd q_old(7); //temporary store of the computed pose
    Eigen::VectorXd nullspace(7); //temporary store of the null space motion part
    Eigen::VectorXd q_diff(7);
    Eigen::VectorXd q(7);
    Eigen::VectorXd q_comf(18); //both arms
    Eigen::VectorXd q_merged(18); //both arms
    float epsilon = 0.001;

    std::cout << "Vars init" << std::endl;

    //Load data & precompute info that is constant in the loop
    bax.GetTargets(target); //load targets (8*3=24 length vector)
    Eigen::VectorXd y = bax.GetIK(qstart1);
    J = bax.GetJ(qstart1);
    Eigen::MatrixXd J_right = J.block(0,0,3,7);
    Eigen::MatrixXd J_left = J.block(6,7,3,7);
    Eigen::MatrixXd Jpinv_right = Winv * J_right.transpose() * (J_right * Winv * J_right.transpose() + Cinv).inverse();
    Eigen::MatrixXd Jpinv_left = Winv * J_left.transpose() * (J_right * Winv * J_left.transpose() + Cinv).inverse();


    /*

    std::cout << "Precomputed info before for" << std::endl;

    //Compute each pose and cost
    for (int i=0; i<8; i++){
        for (int j=0; j<4; j++){
            
            q_comf = (j < 2) ? q_comf1 : q_comf2; //indicate which comfort position to try in this iteration

            if (j%2 == 0){
              //right
              std::cout << "Right arm" << std::endl;
              q = Eigen::VectorXd(qstart1.segment(0,7));
              std::cout << "q: " << q << std::endl;
              while ((q - q_old).cwiseAbs().maxCoeff() > epsilon){
                //Repeat until change is small enough
                q_merged << q, Eigen::VectorXd::Zero(11);
                std::cout << "q_merged: " << q_merged << std::endl;
                bax.SetJointAngles(q_merged);
                bax.AdvanceSimulation();

                y = bax.GetIK(q_merged);
                
                J = bax.GetJ(q_merged);
                J_right = J.block(0,0,3,7);

                std::cout << "Dimensions: " << std::endl;
                std::cout << "Winv " << Winv.rows() << "x" << Winv.cols() << std::endl;
                std::cout << "J_right " << J_right.rows() << "x" << J_right.cols() << std::endl;
                std::cout << "Cinv " << Cinv.rows() << "x" << Cinv.cols() << std::endl;

                Jpinv_right = Winv * J_right.transpose() * (J_right * Winv * J_right.transpose() + Cinv).inverse();
                
                std::cout << "Jpinv_right " << Jpinv_right.rows() << "x" << Jpinv_right.cols() << std::endl;

                nullspace = (Eigen::MatrixXd::Identity(7, 7) - (Jpinv_right * J_right)) * (q_comf.segment(0, 7) - q);  
                
                std::cout << "nullspace " << nullspace.rows() << "x" << nullspace.cols() << std::endl;

                q_diff = (Jpinv_right * (target.segment(i*3,3) - y.segment(0, 3))) + nullspace;
                q_old = q;
                q = q + 0.1 * q_diff;
                //sleep(1);
              }

              q_merged << Eigen::VectorXd::Zero(9), q, Eigen::VectorXd::Zero(2);
              bax.SetJointAngles(q_merged);
              bax.AdvanceSimulation();

              q_diff = qstart1.segment(0,7) - q;

              std::cout << "SUCCESS " << std::endl;
            } else {
              //left
              q = Eigen::VectorXd(qstart1.segment(9,7));
              while ((q - q_old).cwiseAbs().maxCoeff() > epsilon){
                //Repeat until change is small enough
                q_merged << Eigen::VectorXd::Zero(9), q, Eigen::VectorXd::Zero(2);

                bax.SetJointAngles(q_merged);
                bax.AdvanceSimulation();

                y = bax.GetIK(q_merged);

                J = bax.GetJ(q_merged);
                J_left = J.block(6,7,3,7);
                Jpinv_left = Winv * J_left.transpose() * (J_left * Winv * J_left.transpose() + Cinv).inverse();

                nullspace = (Eigen::MatrixXd::Identity(7, 7) - Jpinv_left * J_left) * (q_comf.segment(9, 7) - q);  
                q_diff = Jpinv_left * (target.segment(i*3,3) - y.segment(6, 3)) + nullspace;
                q_old = q;
                q = q + 0.1 * q_diff;
              }

              q_merged << Eigen::VectorXd::Zero(9), q, Eigen::VectorXd::Zero(2);
              bax.SetJointAngles(q_merged);
              bax.AdvanceSimulation();

              q_diff = qstart1.segment(9,7) - q;
            }

            //Compute cost
            costs(i,j) = q_diff.transpose() * Winv.inverse() * q_diff;   
            std::cout << "Cost: " << costs(i,j) << std::endl;         
        }
    }

    std::cout << "Costs matrix\n" << costs << std::endl;

    

    //Part B

    std::cout << "Starting Part B" << std::endl;

    Eigen::VectorXd target0 = target.segment(0,3);

    //1 (comf)
    q = Eigen::VectorXd(qstart1.segment(0,7));
    q_old = Eigen::VectorXd::Zero(7);
    while ((q - q_old).cwiseAbs().maxCoeff() > epsilon){
      //Repeat until change is small enough
      q_merged << q, Eigen::VectorXd::Zero(11);

      bax.SetJointAngles(q_merged);
      bax.AdvanceSimulation();

      y = bax.GetIK(q_merged);
      J = bax.GetJ(q_merged);
      J_right = J.block(0,0,3,7);
      Jpinv_right = Winv * J_right.transpose() * (J_right * Winv * J_right.transpose() + Cinv).inverse();
      

      nullspace = (Eigen::MatrixXd::Identity(7, 7) - Jpinv_right * J_right) * (q_comf.segment(0, 7) - q);  
      q_diff = Jpinv_right * (target0 - y.segment(0, 3)) + nullspace;
      q_old = q;
      q = q + 0.1 * q_diff;
    }

    //2 (comf)
    q = Eigen::VectorXd(qstart2.segment(0,7));
    q_old = Eigen::VectorXd::Zero(7);
    while ((q - q_old).cwiseAbs().maxCoeff() > epsilon){
      //Repeat until change is small enough
      q_merged << q, Eigen::VectorXd::Zero(11);

      bax.SetJointAngles(q_merged);
      bax.AdvanceSimulation();

      y = bax.GetIK(q_merged);
      J = bax.GetJ(q_merged);
      J_right = J.block(0,0,3,7);
      Jpinv_right = Winv * J_right.transpose() * (J_right * Winv * J_right.transpose() + Cinv).inverse();
      

      nullspace = (Eigen::MatrixXd::Identity(7, 7) - Jpinv_right * J_right) * (q_comf.segment(0, 7) - q);  
      q_diff = Jpinv_right * (target0 - y.segment(0, 3)) + nullspace;
      q_old = q;
      q = q + 0.1 * q_diff;
    }

    //3 (comf)
    q = Eigen::VectorXd(qstart3.segment(0,7));
    q_old = Eigen::VectorXd::Zero(7);
    while ((q - q_old).cwiseAbs().maxCoeff() > epsilon){
      //Repeat until change is small enough
      q_merged << q, Eigen::VectorXd::Zero(11);

      bax.SetJointAngles(q_merged);
      bax.AdvanceSimulation();

      y = bax.GetIK(q_merged);
      J = bax.GetJ(q_merged);
      J_right = J.block(0,0,3,7);
      Jpinv_right = Winv * J_right.transpose() * (J_right * Winv * J_right.transpose() + Cinv).inverse();
      

      nullspace = (Eigen::MatrixXd::Identity(7, 7) - Jpinv_right * J_right) * (q_comf.segment(0, 7) - q);  
      q_diff = Jpinv_right * (target0 - y.segment(0, 3)) + nullspace;
      q_old = q;
      q = q + 0.1 * q_diff;
    }

    //1 (no comf)
    q = Eigen::VectorXd(qstart1.segment(0,7));
    q_old = Eigen::VectorXd::Zero(7);
    while ((q - q_old).cwiseAbs().maxCoeff() > epsilon){
      //Repeat until change is small enough
      q_merged << q, Eigen::VectorXd::Zero(11);

      bax.SetJointAngles(q_merged);
      bax.AdvanceSimulation();

      y = bax.GetIK(q_merged);
      J = bax.GetJ(q_merged);
      J_right = J.block(0,0,3,7);
      Jpinv_right = Winv * J_right.transpose() * (J_right * Winv * J_right.transpose() + Cinv).inverse();
      

      //nullspace = (Eigen::MatrixXd::Identity(7, 7) - Jpinv_right * J_right) * (q_comf.segment(0, 7) - q);  
      q_diff = Jpinv_right * (target0 - y.segment(0, 3));// + nullspace;
      q_old = q;
      q = q + 0.1 * q_diff;
    }

    //2 (no comf)
    q = Eigen::VectorXd(qstart2.segment(0,7));
    q_old = Eigen::VectorXd::Zero(7);
    while ((q - q_old).cwiseAbs().maxCoeff() > epsilon){
      //Repeat until change is small enough
      q_merged << q, Eigen::VectorXd::Zero(11);

      bax.SetJointAngles(q_merged);
      bax.AdvanceSimulation();

      y = bax.GetIK(q_merged);
      J = bax.GetJ(q_merged);
      J_right = J.block(0,0,3,7);
      Jpinv_right = Winv * J_right.transpose() * (J_right * Winv * J_right.transpose() + Cinv).inverse();
      

      //nullspace = (Eigen::MatrixXd::Identity(7, 7) - Jpinv_right * J_right) * (q_comf.segment(0, 7) - q);  
      q_diff = Jpinv_right * (target0 - y.segment(0, 3));// + nullspace;
      q_old = q;
      q = q + 0.1 * q_diff;
    }

    //3 (no comf)
    q = Eigen::VectorXd(qstart3.segment(0,7));
    q_old = Eigen::VectorXd::Zero(7);
    while ((q - q_old).cwiseAbs().maxCoeff() > epsilon){
      //Repeat until change is small enough
      q_merged << q, Eigen::VectorXd::Zero(11);

      bax.SetJointAngles(q_merged);
      bax.AdvanceSimulation();

      y = bax.GetIK(q_merged);
      J = bax.GetJ(q_merged);
      J_right = J.block(0,0,3,7);
      Jpinv_right = Winv * J_right.transpose() * (J_right * Winv * J_right.transpose() + Cinv).inverse();
      

      //nullspace = (Eigen::MatrixXd::Identity(7, 7) - Jpinv_right * J_right) * (q_comf.segment(0, 7) - q);  
      q_diff = Jpinv_right * (target0 - y.segment(0, 3));// + nullspace;
      q_old = q;
      q = q + 0.1 * q_diff;
    }

    */


    std::cout << "PART C" << std::endl;
    sleep(5);

    // PART C
    Eigen::MatrixXd q_mat(24,7);

    for (int i=0; i<8; i++){
      for (int j=0; j<3; j++){
        std::cout << "New iteration, i: " << i << "   j: " << j << std::endl;
        switch (j){
          case 0:
              q = Eigen::VectorXd(qstart1.segment(0,7));
            break;
          case 1:
              q = Eigen::VectorXd(qstart2.segment(0,7));
            break;
          case 2:
              q = Eigen::VectorXd(qstart3.segment(0,7));
            break;
        };  
        q_old = Eigen::VectorXd::Zero(7);
        while ((q - q_old).cwiseAbs().maxCoeff() > epsilon){
          //Repeat until change is small enough
          q_merged << q, Eigen::VectorXd::Zero(11);
          //std::cout << "q_merged: " << q_merged << std::endl;
          bax.SetJointAngles(q_merged);
          bax.AdvanceSimulation();

          y = bax.GetIK(q_merged);
          
          J = bax.GetJ(q_merged);
          J_right = J.block(0,0,3,7);

          /*std::cout << "Dimensions: " << std::endl;
          std::cout << "Winv " << Winv.rows() << "x" << Winv.cols() << std::endl;
          std::cout << "J_right " << J_right.rows() << "x" << J_right.cols() << std::endl;
          std::cout << "Cinv " << Cinv.rows() << "x" << Cinv.cols() << std::endl;
          */
          
          Jpinv_right = Winv * J_right.transpose() * (J_right * Winv * J_right.transpose() + Cinv).inverse();
          
          //std::cout << "Jpinv_right " << Jpinv_right.rows() << "x" << Jpinv_right.cols() << std::endl;

          //nullspace = (Eigen::MatrixXd::Identity(7, 7) - (Jpinv_right * J_right)) * (q_comf.segment(0, 7) - q);  
          
          //std::cout << "nullspace " << nullspace.rows() << "x" << nullspace.cols() << std::endl;

          q_diff = (Jpinv_right * (target.segment(i*3,3) - y.segment(0, 3))); //+ nullspace;
          q_old = q;
          q = q + 0.1 * q_diff;
          //sleep(1);
        }
        std::cout << "Assigning row to q_mat" << std::endl;
        q_mat.block(i*3+j,0,1,7) = q.transpose();
        std::cout << "Assigned row to q_mat" << std::endl;
      }
    }

    //TO DO: RUN PCA
    Eigen::MatrixXd centered = q_mat.rowwise() - q_mat.colwise().mean();
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(centered, Eigen::ComputeThinU | Eigen::ComputeThinV);
    //JacobiSVD<MatrixXf> svd(centered, ComputeFullU | ComputeFullV);
    std::cout << "Singular values:" << std::endl << svd.singularValues() << std::endl;
    std::cout << "Singular vectors:" << std::endl << svd.matrixU() << std::endl;
    
    sleep(20);

    //DEMO
    std::cout << "Right arm" << std::endl;
    q = Eigen::VectorXd(qstart1.segment(0,7));
    std::cout << "q: " << q << std::endl;
    q_old = Eigen::VectorXd::Zero(7);
    q_comf = q_comf1;
    while ((q - q_old).cwiseAbs().maxCoeff() > epsilon){
      //Repeat until change is small enough
      q_merged << q, Eigen::VectorXd::Zero(11);
      std::cout << "q_merged: " << q_merged << std::endl;
      bax.SetJointAngles(q_merged);
      bax.AdvanceSimulation();

      y = bax.GetIK(q_merged);
      
      J = bax.GetJ(q_merged);
      J_right = J.block(0,0,3,7);

      std::cout << "Dimensions: " << std::endl;
      std::cout << "Winv " << Winv.rows() << "x" << Winv.cols() << std::endl;
      std::cout << "J_right " << J_right.rows() << "x" << J_right.cols() << std::endl;
      std::cout << "Cinv " << Cinv.rows() << "x" << Cinv.cols() << std::endl;

      Jpinv_right = Winv * J_right.transpose() * (J_right * Winv * J_right.transpose() + Cinv).inverse();
      
      std::cout << "Jpinv_right " << Jpinv_right.rows() << "x" << Jpinv_right.cols() << std::endl;

      nullspace = (Eigen::MatrixXd::Identity(7, 7) - (Jpinv_right * J_right)) * (q_comf.segment(0, 7) - q);  
      
      std::cout << "nullspace " << nullspace.rows() << "x" << nullspace.cols() << std::endl;

      q_diff = (Jpinv_right * (target.segment(7*3,3) - y.segment(0, 3))) + nullspace;
      q_old = q;
      q = q + 0.1 * q_diff;
      //sleep(1);
    }    


    ////////////////////////////////////////////////

    // Update simulation
    bax.AdvanceSimulation();
  }
  // Stop simulation and close connection
  bax.StopSimulation();
	return(0);
}

