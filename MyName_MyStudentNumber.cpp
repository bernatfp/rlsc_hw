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
                
                J = bax.GetJ(qstart1);
                J_right = J.block(0,0,3,7);

                std::cout << "Dimensions: " << std::endl;
                std::cout << "Winv " << Winv.size() << std::endl;
                std::cout << "J_right " << J_right.size() << std::endl;
                std::cout << "Cinv " << Cinv.size() << std::endl;

                Jpinv_right = Winv * J_right.transpose() * (J_right * Winv * J_right.transpose() + Cinv).inverse();
                
                std::cout << "Jpinv_right " << Jpinv_right.size() << std::endl;

                nullspace = (Eigen::MatrixXd::Identity(J_right.rows(), J_right.cols()) - Jpinv_right * J_right) * (q_comf.segment(0, 7) - q);  
                
                std::cout << "nullspace " << nullspace.size() << std::endl;

                q_diff = Jpinv_right * (target.segment(i*3,3) - y.segment(0, 3)) + nullspace;
                q_old = q;
                q = q + q_diff;
              }
            } else {
              //left
              q = Eigen::VectorXd(qstart1.segment(9,7));
              while ((q - q_old).cwiseAbs().maxCoeff() > epsilon){
                //Repeat until change is small enough
                q_merged << Eigen::VectorXd::Zero(9), q, Eigen::VectorXd::Zero(2);

                bax.SetJointAngles(q_merged);
                bax.AdvanceSimulation();

                y = bax.GetIK(q_merged);

                J = bax.GetJ(qstart1);
                J_left = J.block(6,7,3,7);
                Jpinv_left = Winv * J_left.transpose() * (J_right * Winv * J_left.transpose() + Cinv).inverse();

                nullspace = (Eigen::MatrixXd::Identity(J_right.rows(), J_right.cols()) - Jpinv_left * J_left) * (q_comf.segment(9, 7) - q);  
                q_diff = Jpinv_left * (target.segment(i*3,3) - y.segment(6, 3)) + nullspace;
                q_old = q;
                q = q + q_diff;
              }
            }

            //Compute cost
            q_diff = q - qstart1;
            costs(i,j) = q_diff.transpose() * Winv.inverse() * q_diff;            
        }
    }

    std::cout << "Costs matrix\n" << costs << std::endl;


    //Part B

    Eigen::VectorXd target0 = target.segment(0,3);

    //1 (comf)
    q = Eigen::VectorXd(qstart1.segment(0,7));
    q_old = Eigen::VectorXd::Zero(7);
    while ((q - q_old).cwiseAbs().maxCoeff() > epsilon){
      //Repeat until change is small enough
      q_merged << q, Eigen::VectorXd::Zero(11);
      y = bax.GetIK(q_merged);
      
      J = bax.GetJ(qstart1);
      J_right = J.block(0,0,3,7);
      Jpinv_right = Winv * J_right.transpose() * (J_right * Winv * J_right.transpose() + Cinv).inverse();
      

      nullspace = (Eigen::MatrixXd::Identity(J_right.rows(), J_right.cols()) - Jpinv_right * J_right) * (q_comf.segment(0, 7) - q);  
      q_diff = Jpinv_right * (target0 - y.segment(0, 3)) + nullspace;
      q_old = q;
      q = q + q_diff;
    }

    //2 (comf)
    q = Eigen::VectorXd(qstart2.segment(0,7));
    q_old = Eigen::VectorXd::Zero(7);
    while ((q - q_old).cwiseAbs().maxCoeff() > epsilon){
      //Repeat until change is small enough
      q_merged << q, Eigen::VectorXd::Zero(11);
      y = bax.GetIK(q_merged);
      
      J = bax.GetJ(qstart1);
      J_right = J.block(0,0,3,7);
      Jpinv_right = Winv * J_right.transpose() * (J_right * Winv * J_right.transpose() + Cinv).inverse();
      
      nullspace = (Eigen::MatrixXd::Identity(J_right.rows(), J_right.cols()) - Jpinv_right * J_right) * (q_comf.segment(0, 7) - q);  
      q_diff = Jpinv_right * (target0 - y.segment(0, 3)) + nullspace;
      q_old = q;
      q = q + q_diff;
    }

    //3 (comf)
    q = Eigen::VectorXd(qstart3.segment(0,7));
    q_old = Eigen::VectorXd::Zero(7);
    while ((q - q_old).cwiseAbs().maxCoeff() > epsilon){
      //Repeat until change is small enough
      q_merged << q, Eigen::VectorXd::Zero(11);
      y = bax.GetIK(q_merged);
      
      J = bax.GetJ(qstart1);
      J_right = J.block(0,0,3,7);
      Jpinv_right = Winv * J_right.transpose() * (J_right * Winv * J_right.transpose() + Cinv).inverse();
      
      nullspace = (Eigen::MatrixXd::Identity(J_right.rows(), J_right.cols()) - Jpinv_right * J_right) * (q_comf.segment(0, 7) - q);  
      q_diff = Jpinv_right * (target0 - y.segment(0, 3)) + nullspace;
      q_old = q;
      q = q + q_diff;
    }

    //1 (no comf)
    q = Eigen::VectorXd(qstart1.segment(0,7));
    q_old = Eigen::VectorXd::Zero(7);
    while ((q - q_old).cwiseAbs().maxCoeff() > epsilon){
      //Repeat until change is small enough
      q_merged << q, Eigen::VectorXd::Zero(11);
      y = bax.GetIK(q_merged);
      
      J = bax.GetJ(qstart1);
      J_right = J.block(0,0,3,7);
      Jpinv_right = Winv * J_right.transpose() * (J_right * Winv * J_right.transpose() + Cinv).inverse();
      
      //nullspace = (Eigen::MatrixXd::Identity(J_right.size()) - Jpinv_right * J_right) * (q_comf.segment(0, 7) - q);  
      q_diff = Jpinv_right * (target0 - y.segment(0, 3)); // + nullspace;
      q_old = q;
      q = q + q_diff;
    }

    //2 (no comf)
    q = Eigen::VectorXd(qstart2.segment(0,7));
    q_old = Eigen::VectorXd::Zero(7);
    while ((q - q_old).cwiseAbs().maxCoeff() > epsilon){
      //Repeat until change is small enough
      q_merged << q, Eigen::VectorXd::Zero(11);
      y = bax.GetIK(q_merged);
      
      J = bax.GetJ(qstart1);
      J_right = J.block(0,0,3,7);
      Jpinv_right = Winv * J_right.transpose() * (J_right * Winv * J_right.transpose() + Cinv).inverse();
      
      //nullspace = (Eigen::MatrixXd::Identity(J_right.size()) - Jpinv_right * J_right) * (q_comf.segment(0, 7) - q);  
      q_diff = Jpinv_right * (target0 - y.segment(0, 3)); // + nullspace;
      q_old = q;
      q = q + q_diff;
    }

    //3 (no comf)
    q = Eigen::VectorXd(qstart3.segment(0,7));
    q_old = Eigen::VectorXd::Zero(7);
    while ((q - q_old).cwiseAbs().maxCoeff() > epsilon){
      //Repeat until change is small enough
      q_merged << q, Eigen::VectorXd::Zero(11);
      y = bax.GetIK(q_merged);
      
      J = bax.GetJ(qstart1);
      J_right = J.block(0,0,3,7);
      Jpinv_right = Winv * J_right.transpose() * (J_right * Winv * J_right.transpose() + Cinv).inverse();
      
      //nullspace = (Eigen::MatrixXd::Identity(J_right.size()) - Jpinv_right * J_right) * (q_comf.segment(0, 7) - q);  
      q_diff = Jpinv_right * (target0 - y.segment(0, 3)); // + nullspace;
      q_old = q;
      q = q + q_diff;
    }






    // PART C



    ////////////////////////////////////////////////

    // Update simulation
    bax.AdvanceSimulation();
  }
  // Stop simulation and close connection
  bax.StopSimulation();
	return(0);
}

