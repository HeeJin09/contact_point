#include <stdio.h>
#include <signal.h>
#include <unistd.h>
#include <sys/mman.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <native/task.h>
#include <native/timer.h>
#include <rtdk.h>
#include <../eigen/Eigen/Dense>
#include "NRMKsercan_tp.h"
#include "NRMKhw_tp.h"
#include <vector>
#include <cmath>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <netdb.h>
#include <errno.h>
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <unistd.h>
#include <arpa/inet.h>
#include <sys/socket.h>
#include <fstream>  
#include <Eigen/Dense>
#define PORT 8000

using namespace std;
using namespace Eigen;


TP_PORT test_port=1;

int rtsercan_fd  = -1;

RT_TASK write_task;
RT_TASK read_task;
int working=1;

//task start value
double trigger = 0;
int fd=-1;
void * read_thread(void* arg);
//CANbus data to physical data
unsigned char data_field[16]; // storage buffer for data field
//..... Received CAN data Save .....
// 8 byte data of can message id is #1 save in data_field [0] ~ [7]
// 8 byte data of can message id is #2 save in data_field [8] ~ [15] // data field processing
short raw_data[6] = { 0 };
unsigned short temp;
unsigned DF=50, DT=2000; // DF, DT depend on the model, refer to 3.6.11
double ft_array[6];
double ft_array_init[6];

double fx = 0;
double fy = 0;
double fz = 0;
double mx = 0;
double my = 0;
double mz = 0;
double norm_f = 0;
double ux = 0;
double uy = 0;
double uz = 0;
double filtered_fx = 0.0;
double filtered_fy = 0.0;
double filtered_fz = 0.0;
double filtered_mx = 0.0;
double filtered_my = 0.0;
double filtered_mz = 0.0;

double filter_f[6];
double x[3];
double pre_filter_f[6] ={0,0,0,0,0,0};
double g_force[6];
bool init_flag = 0;

void touch_point(double* x, int size);
using namespace Eigen;
using namespace std;
double qx = 0.0, qy = 0.0, qz = 0.0, qw = 0.0; 
Eigen::Matrix3d R_0;
int count_Test = 0;
int count_Test_1 = 0;
double data_1[6][1000];
int sock = 0;
struct sockaddr_in serv_addr; 
char *ip = "192.168.0.19";
tuple<double, double, double, double, double, double> moving_average_filter(
    double fx, double fy, double fz, double mx, double my, double mz)
{
    int window_size = 100;
    static vector<double> input_fx;
    static vector<double> input_fy;
    static vector<double> input_fz;
    static vector<double> input_mx;
    static vector<double> input_my;
    static vector<double> input_mz;

    input_fx.push_back(fx);
    input_fy.push_back(fy);
    input_fz.push_back(fz);
    input_mx.push_back(mx);
    input_my.push_back(my);
    input_mz.push_back(mz);

    if (input_fx.size() >= window_size)
    {
        double sum_fx = 0;
        double sum_fy = 0;
        double sum_fz = 0;
        double sum_mx = 0;
        double sum_my = 0;
        double sum_mz = 0;

        for (int i = input_fx.size() - window_size; i < input_fx.size(); i++)
        {
            sum_fx += input_fx[i];
            sum_fy += input_fy[i];
            sum_fz += input_fz[i];
            sum_mx += input_mx[i];
            sum_my += input_my[i];
            sum_mz += input_mz[i];
        }

        double filtered_value_fx = sum_fx / window_size;
        double filtered_value_fy = sum_fy / window_size;
        double filtered_value_fz = sum_fz / window_size;
        double filtered_value_mx = sum_mx / window_size;
        double filtered_value_my = sum_my / window_size;
        double filtered_value_mz = sum_mz / window_size;

        input_fx.erase(input_fx.begin());
        input_fy.erase(input_fy.begin());
        input_fz.erase(input_fz.begin());
        input_mx.erase(input_mx.begin());
        input_my.erase(input_my.begin());
        input_mz.erase(input_mz.begin());

        return {filtered_value_fx, filtered_value_fy, filtered_value_fz,
                filtered_value_mx, filtered_value_my, filtered_value_mz};
    }

    return {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
}
double filter_f__[6] = {0.0,0.0,0.0,0.0,0.0,0.0};
double pre_filter_f__[6] ={0.0,0.0,0.0,0.0,0.0,0.0};

tuple<double, double, double, double, double, double> LPF_filter(
    double fx, double fy, double fz, double mx, double my, double mz)
{
   double Ts= 1/DT;
    VectorXd x_k(6); 
    x_k << fx, fy, fz, mx, my, mz ;// Measurement vector

    for (int n = 0; n < 6; n++)
    {
            float alpha = 0.909;
            filter_f__[n] = alpha * pre_filter_f__[n] + (1-alpha)*(x_k[n]);
            pre_filter_f__[n] = filter_f__[n];
    } 
    return {filter_f__[0], filter_f__[1], filter_f__[2], filter_f__[3], filter_f__[4], filter_f__[5]};
}
Eigen::Matrix<double, 6, 6> Pk_prev = Matrix<double, 6, 6>::Identity();
Eigen::Matrix<double, 6, 1> xk_prev = Matrix<double, 6, 1>::Zero();

tuple<double, double, double, double, double, double> Kalman_filter(
    double fx, double fy, double fz, double mx, double my, double mz)
{
   
    VectorXd x(6); // State vector
    MatrixXd P(6, 6); // Error covariance matrix
    MatrixXd F(6, 6); // State transition matrix
    MatrixXd Q(6, 6); // Process noise covariance matrix
    MatrixXd H(6, 6); // Measurement matrix
    MatrixXd R(6, 6); 
    R << 0.09, 0, 0, 0, 0, 0,
    0, 0.04, 0, 0, 0, 0,
    0, 0, 0.16, 0, 0, 0, 
    0, 0, 0, 0.001, 0, 0,
    0, 0, 0, 0, 0.001, 0,
    0, 0, 0, 0, 0, 0.001;// Measurement noise covariance matrix
    VectorXd y(6); 
    y << fx, fy, fz, mx, my, mz ;// Measurement vector
    VectorXd x_hat_minus(6); // Predicted state vector
    MatrixXd Pk_minus(6, 6); // Predicted error covariance matrix
    MatrixXd K(6, 6); // Kalman gain matrix
    VectorXd x_hat_plus(6); // Updated state vector
    MatrixXd Pk_plus(6, 6); // Updated error covariance matrix

    // Set initial values
    x.setZero(); // Set initial state to zero
    P.setIdentity(); // Set initial error covariance to identity matrix
    F.setIdentity(); // Set state transition matrix to identity matrix
    Q.setIdentity(); // Set process noise covariance to identity matrix
    H.setIdentity(); // Set measurement matrix to identity matrix
    R.setIdentity(); // Set measurement noise covariance to identity matrix
    MatrixXd teststests(6, 6); 
    teststests << 0.3, 0, 0, 0, 0, 0,
    0, 0.3, 0, 0, 0, 0,
    0, 0, 0.3, 0, 0, 0, 
    0, 0, 0, 0.3, 0, 0,
    0, 0, 0, 0, 0.3, 0,
    0, 0, 0, 0, 0, 0.3;
    // Perform Kalman filter

        // Prediction step
    x_hat_minus = F * x;
    Pk_minus = F * P * F.transpose() + Q;

    // Update step
    K = Pk_minus * H.transpose() * (H * Pk_minus * H.transpose() + R).inverse() ;
    x_hat_plus = x_hat_minus + K * (y - H * x_hat_minus);
    Pk_plus = (MatrixXd::Identity(6, 6) - K * H) * Pk_minus * (MatrixXd::Identity(6, 6) - K * H).transpose() + K * R * K.transpose();

    // Update state and error covariance
    x = x_hat_plus;
    P = Pk_plus;

    // Print results

   //cout << "K  = " << K .transpose() << std::endl;
   //cout << "Pk_plus = " << Pk_plus << std::endl;

    return {x_hat_plus[0], x_hat_plus[1], x_hat_plus[2], x_hat_plus[3], x_hat_plus[4], x_hat_plus[5]};
    
}
void touch_point(double* x, int size)
{
    int t_rx, t_ry;
    if(int(x[1]*100) < 0) t_rx = (int)(x[1] *100) + 6;
    else if (int(abs(x[1]*100))== 0) t_rx = 6;
    else if (int(abs(x[1]*100)) > 0)  t_rx = (int)(x[1] *100) + 6;

    if(int(x[0]*100) < 0) t_ry = abs((int)(x[0] *100)) + 6;
    else if (int(abs(x[0]*100))== 0) t_ry =6 ;
    else if (int(abs(x[0]*100)) > 0)  t_ry = 6 - abs((int)(x[0] *100));

    for(int i = 0; i < 11 ; i++)
    {
        for(int j = 0; j < 11; j++)
        {
            if(abs(t_rx) == (i+1) && abs(t_ry) == (j+1)) 
            {
                printf("  ");
            }
            else  
            {
                printf("██");
            }
        }
        printf("\n");
    }

    printf("\n");
    printf("\n");
}
void * read_thread(void* arg)
{
    int nbytes;
    char chr[100];
    vector<double> vec(4);
    while(1)
    {
        
        vector<double> input_qx;
        vector<double> input_qy;
        vector<double> input_qz;
        vector<double> input_qw;

        int window_size = 100;
        nbytes=read(fd,&chr,100);
        if (nbytes>0){
            // ','로 구분된 문자열을 토큰으로 분리하여 Eigen Vector에 저장
            double filtered_value_qx = 0;
            double filtered_value_qy = 0;
            double filtered_value_qz = 0;
            double filtered_value_qw = 0;
            char *token = strtok(chr, ",");
            for (int i = 1; i <= 5 && token != NULL; i++) {
                if (i >= 2 && i <= 5) {
                    vec[i-2] = atof(token);
                }
                token = strtok(NULL, ",");
            }          
            qx = vec[0];
            qy = vec[1];
            qz = vec[2];
            qw = vec[3];

            Eigen::Quaterniond q_1(qx, qy, qz, qw);

            Eigen::Matrix3d R_;
                R_ <<  -1.0, 0.0, 0.0,
                    0.0, 1.0, 0.0,
                    0.0, 0.0 , 1.0;
            Eigen::Matrix3d rotation_matrix = q_1.toRotationMatrix();  // 쿼터니언을 회전 행렬로 변환합니다.
            Eigen::Matrix3d R_trans = R_.transpose();
            Eigen::Matrix3d result_rotation_matrix = R_ * rotation_matrix;
            if(count_Test == 0){
                    R_0 =  result_rotation_matrix;      
                    count_Test = 1;
                }
        }
        usleep(1e4);
    }
}

void cleanup_all(void)
{
    working=0;
    rt_task_delete(&write_task);
    rt_task_delete(&read_task);

}

void catch_signal(int sig)
//void cath_signal()
{
    cleanup_all();
    return;
}


void write_thread(void* arg)
{
    CAN_FRAME cfame;

    rt_task_set_periodic(NULL, TM_NOW, 1e9);
    while (working)
    {
        //initial setting
        //rt_printf("writing task initial setting completed...  \n");
        cfame.can_id=0x64;
        cfame.can_dlc=8;
        for (int i=0; i<8; ++i) cfame.data[i]=0x00;

        //start sensing
        if(trigger == 0)
        {
            //rt_printf("start sensing...   \n");
            for (int i=0; i<8; ++i){
            if(i==0)
                cfame.data[i]=0x0B;
            else
                cfame.data[i]=0x00;
        }
    }

        RTSERCAN_write(rtsercan_fd, cfame);
        rt_task_wait_period(NULL);
    }

}
Eigen::Matrix3d skew(Eigen::Matrix<double,3,1> w){
    Eigen::Matrix3d ret=Eigen::Matrix3d::Zero();
    ret(0,1) = -w(2);
    ret(0,2) = w(1);
    ret(1,0) = w(2);
    ret(1,2) = -w(0);
    ret(2,0) = -w(1);
    ret(2,1) = w(0);
    return ret;
}
Eigen::Matrix<double,6,6> Ad(Eigen::Matrix3d R, Eigen::Matrix<double,3,1> p){
    Eigen::Matrix<double,6,6> ret;
    ret<<0,0,0,0,0,0,0,0,0,0,0,0,
         0,0,0,0,0,0,0,0,0,0,0,0,
         0,0,0,0,0,0,0,0,0,0,0,0;
    ret.block<3,3>(0,0) = R;
    ret.block<3,3>(3,3) = R;
    ret.block<3,3>(3,0) = skew(p)*R;
    return ret;
}
void sock_Test()
{
    if ((sock = socket(AF_INET, SOCK_DGRAM, 0)) < 0) 
        { 
            printf("\n Socket creation error \n"); 
            return -1; 
        } 

        serv_addr.sin_family = AF_INET; 
        serv_addr.sin_port = htons(PORT); 


        // Convert IPv4 and IPv6 addresses from text to binary form 
        if(inet_pton(AF_INET, "192.168.0.19", &serv_addr.sin_addr)<=0)  
        { 
            printf("\nInvalid address/ Address not supported \n"); 
            return -1; 
        } 

        if (connect(sock, (struct sockaddr *)&serv_addr, sizeof(serv_addr)) < 0) 
        { 
            printf("\nConnection Failed \n"); 
            return -1; 
        } 
}
void read_task_proc(void *arg)
{
    int res1, res2;
    CAN_FRAME RxFrame1;
    CAN_FRAME RxFrame2;
    rt_task_set_periodic(NULL, TM_NOW, 1e6);
    static unsigned int print_count = 0;
    std::ofstream file("data.csv", std::ios::out); // 파일 열기

    if (!file.is_open()) {
        std::cerr << "Failed to open file!" << std::endl;
        return 1;
    }

    while (working)
    {
        //Get CANbus Torque value
        res1=RTSERCAN_read(rtsercan_fd, &RxFrame1);
        res2=RTSERCAN_read(rtsercan_fd, &RxFrame2);


        if (res1==SERCAN_ERR_FREE && res2==SERCAN_ERR_FREE)
        {/*
            rt_print_CANFrame(RxFrame1);
            rt_print_CANFrame(RxFrame2);
        */

        //CANbus data to Torque data
        for(int i=0; i<8; ++i)
        {
            data_field[i] = (unsigned char) RxFrame1.data[i];
            data_field[i+8] = (unsigned char) RxFrame2.data[i];
        }

        for (int idx = 0; idx < 6; idx++)
        {
            temp = data_field [2 * idx + 1] * 256;
            temp += data_field [2 * idx + 2];

            raw_data[idx] = (signed short)temp; // variable casting
        }
        // Conversion from signed short data to float data and data scaling

        // Set Force/Torque Original
        for (int n = 0; n < 3; n++)
        {
            ft_array[n] = (((double)raw_data[n] ) / DF); // refer to 3.6.11
            ft_array[n + 3] = (((double)raw_data[n + 3] ) / DT); // refer to 3.6.11
        }
        int num_x, num_y, num_z;
        if(init_flag == 0){
            for(int i=0;i<6;i++){
                ft_array_init[i] = ft_array[i]; 
            }
            
            init_flag = 1;
        }   
        double fx = ft_array[0], fy = ft_array[1], fz = ft_array[2];
        double mx = ft_array[3], my = ft_array[4], mz = ft_array[5];

        double filtered_fx, filtered_fy, filtered_fz, filtered_mx, filtered_my, filtered_mz;
        double k_iltered_fx, k_filtered_fy, k_filtered_fz, k_filtered_mx, k_filtered_my, k_filtered_mz;
        double l_iltered_fx, l_filtered_fy, l_filtered_fz, l_filtered_mx, l_filtered_my, l_filtered_mz;
        
        if (count_Test_1 < 100)
        {
            tie(filtered_fx, filtered_fy, filtered_fz, filtered_mx, filtered_my, filtered_mz) = moving_average_filter(fx, fy, fz, mx, my, mz);
            
        }


        Eigen::Matrix<double,6,1> Wrench_init;
        Wrench_init<<0.0001,
        0.0001,
        0.0001,
        0.0,
        0.0,
       -1.1;

        Eigen::Matrix<double,6,1> Wrench;
        Wrench<<mx,
        my,
        mz,
        fx,
        fy,
        fz;

        Eigen::Matrix<double,6,1> base;
        base<<filtered_mx,
        filtered_my,
        filtered_mz,
        filtered_fx,
        filtered_fy,
        filtered_fz +1.1;
        Eigen::VectorXd ft_finit;
        
        //double mx_1 = ft_finit[0], my_1 = ft_finit[1], mz_1 = ft_finit[2], fx_1 = ft_finit[3], fy_1 = ft_finit[4], fz_1 = ft_finit[5];
        
        //tie(k_iltered_fx, k_filtered_fy, k_filtered_fz, k_filtered_mx, k_filtered_my, k_filtered_mz) = Kalman_filter(fx_1, fy_1, fz_1, mx_1, my_1, mz_1);
        //tie(l_iltered_fx, l_filtered_fy, l_filtered_fz, l_filtered_mx, l_filtered_my, l_filtered_mz) = LPF_filter(fx_1, fy_1, fz_1, mx_1, my_1, mz_1);
                //////////////////////////////////////////////////////////
        
        if ( count_Test_1 > 101)
        {   
            ft_finit = Wrench -  base;
            double mx_1 = ft_finit[0], my_1 = ft_finit[1], mz_1 = ft_finit[2], fx_1 = ft_finit[3], fy_1 = ft_finit[4], fz_1 = ft_finit[5];
            
            //tie(k_iltered_fx, k_filtered_fy, k_filtered_fz, k_filtered_mx, k_filtered_my, k_filtered_mz) = Kalman_filter(fx_1, fy_1, fz_1, mx_1, my_1, mz_1);
            tie(l_iltered_fx, l_filtered_fy, l_filtered_fz, l_filtered_mx, l_filtered_my, l_filtered_mz) = LPF_filter(ft_finit[0], ft_finit[1], ft_finit[2], ft_finit[3], ft_finit[4],  ft_finit[5]);
            
            //data_1[0][count_Test_1 - 101] = l_filtered_mx;
            //data_1[1][count_Test_1 - 101] = l_filtered_my;
            //data_1[2][count_Test_1 - 101] = l_filtered_mz;
            //data_1[3][count_Test_1 - 101] = l_iltered_fx;
            //data_1[4][count_Test_1 - 101] = l_filtered_fy;
            //data_1[5][count_Test_1 - 101] = l_filtered_fz;
            
        }
        //if (count_Test_1 == 1101) {
        //    cout << "start save" << endl;
        //    for (int i = 0; i < 1000; ++i) {
        //        file << data_1[0][i] << ',' << data_1[1][i] << ',' << data_1[2][i] << ','
        //             << data_1[3][i] << ',' << data_1[4][i] << ',' << data_1[5][i] << '\n';
        //    }
        //    file.close(); // 파일 닫기
        //    cout << "file closed" << endl;
        //    break;
        //}

       /////////////////////////////////////////////////////
        count_Test_1 = count_Test_1+1;
        
        
        //Eigen::Matrix<double,6,1> kalman;
        //kalman<<k_iltered_fx,
        //k_filtered_fy,
        //k_filtered_fz,
        //k_filtered_mx,
        //k_filtered_my,
        //k_filtered_mz;
        
        Eigen::Matrix<double,6,1> LPF;
        LPF<<
        l_filtered_mx,
        l_filtered_my,
        l_filtered_mz,
        l_iltered_fx,
        l_filtered_fy,
        l_filtered_fz;



        Eigen::Vector3d f(LPF[3], LPF[4], LPF[5]);
        Eigen::Vector3d m(LPF[0], LPF[1], LPF[2]);
        Eigen::Quaterniond q(qx, qy, qz, qw);

        Eigen::Matrix3d R_;
            R_ << -1.0, 0.0, 0.0,
                    0.0, 1.0, 0.0,
                    0.0, 0.0 , 1.0;
        Eigen::Matrix3d rotation_matrix = q.toRotationMatrix();  // 쿼터니언을 회전 행렬로 변환합니다.
        Eigen::Matrix3d result_rotation_matrix = R_* rotation_matrix;
        Eigen::Matrix3d rotation_matrix_1 = R_0.transpose() * result_rotation_matrix ; 
        
        Eigen::Vector3d p_b(0,0,0.01);
        Eigen::Vector3d p_w(0,0,-0.01);

        
        
        Eigen::VectorXd adjoint_force =  (Ad(rotation_matrix_1,p_b).transpose() * Wrench_init);
        Eigen::Vector3d force_Test (adjoint_force[3], adjoint_force[4], adjoint_force[5]);
        Eigen::Vector3d moment_test = skew(p_b) * force_Test; 
        Eigen::Matrix<double,6,1> body_force;
        body_force<<moment_test[0],
            moment_test[1],
            moment_test[2],
            adjoint_force[3],
            adjoint_force[4],
            adjoint_force[5] + 1.1;
        Eigen::VectorXd body_frame_force = LPF - body_force;

    
        double nx = 0;
        double ny = 0;
        double nz = 1;
        Eigen::Vector3d n(nx, ny, nz);

        double d = 0.01;

        norm_f = f.norm();

        ux = (f[0]/norm_f);
        uy = (f[1]/norm_f);
        uz = (f[2]/norm_f);

        Eigen::Matrix<double,4,3> A;
        A<<0, -uz, uy,
            uz, 0, -ux,
            -uy, ux, 0,
            nx, ny, nz;
        
        Eigen::Matrix<double,4,1> b;
        b<<-1*m[0]/norm_f,
            -1*(m[1]/norm_f),
            -1*(m[2]/norm_f),
             d;

        
        //Eigen::Vector3d world_frame_moment = (rotation_matrix_1.transpose() * body_frame_moment);
        Eigen::MatrixXd A_transpose = A.transpose();
        Eigen::MatrixXd A_transpose_A = A_transpose * A;
        Eigen::MatrixXd A_transpose_A_inv = A_transpose_A.inverse();
        Eigen::VectorXd x = A_transpose_A_inv * A_transpose * b; 
       
        struct Data {
            double matrix[3][3];
            double f[3];
            double m[3];
        };
        Data my_data;
        

        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                my_data.matrix[i][j] = rotation_matrix_1(i,j);
            }
            my_data.m[i] = body_frame_force[i];
            my_data.f[i] = body_frame_force[i+3];
        }
       if (count_Test_1 > 101)
         {   if(print_count % 50 == 0){        
                system("clear");
                cout.precision(6);
                sock_Test();
                send(sock, (char *)&my_data, sizeof(my_data), 0);
                cout << "x = " << x(0) <<",\t" <<x(1) << ",\t" <<x(2) << std::endl;
                cout << " body force: " <<  base.transpose() << endl;
                cout << " adjoint force: " << body_force.transpose() << endl;
                cout << " error force: " <<  body_frame_force.transpose() << endl;
                cout << "rotation matirx  : \n" << rotation_matrix_1 << endl;
                if(norm_f > 5) touch_point(x.data(), x.size()); 
            }
            
        }
            print_count++;
        }
        rt_task_wait_period(NULL);
    }
    //file.close();

}

int main(int argc, char* argv[])
{    
    signal(SIGTERM, catch_signal);
    signal(SIGINT, catch_signal);

    /* no memory-swapping for this programm */
    mlockall(MCL_CURRENT | MCL_FUTURE);
    printf("preventing memory-swapping for this program...  \n");

    // Perform auto-init of rt_print buffers if the task doesn't do so
    rt_print_auto_init(1);

    // open rtsercan*******************
    rtsercan_fd=RTSERCAN_open();
    if (rtsercan_fd < 0) return 1;

    printf("rtsercan is opened! \n");
    //********************************
   
    //-------------------------------------------------------

    printf("creating RT_task... \n");

    rt_task_create(&write_task, "write_task", 0, 99, 0);
    rt_task_start(&write_task, &write_thread, NULL);

    printf("write task created... \n");

    rt_task_create(&read_task, "read_task", 0, 99, 0);
    rt_task_start(&read_task, &read_task_proc, NULL);
    

    cout << "test " << std::endl;
    printf("read task created... \n");
    fd=tp_open_serial_port("/dev/ttyUSB0", 115200);
	pthread_t tid1,tid2;
	pthread_create(&tid1, NULL, read_thread, NULL);
    close(sock);

    pause();

return 0;

}
