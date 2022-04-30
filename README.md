# NaveCodePro
INS/GNSS Integration navigation; Contain Inpure navigation， Integration navigation，Transfer Alignment，Robuster Filter,NHC,etc.

# NaveCodePro使用说明

1，本代码工具箱需要在严恭敏老师的PSINS工具箱的基础上进行使用

​		PSINS相关资料参考严老师网站：http://www.psins.org.cn/

2，NaveCodePro在使用时需要所有文件全部选中在当前路径，以防缺少函数

3，NaveCode后期的维护和更新均上传至GitHub：https://github.com/zelanzou/NaveCodePro

​	  部分数据由于过大，上传至百度网盘：链接：https://pan.baidu.com/s/1VdDXrftKy3uswVjWyirzbQ ，提取码：oajl

# NaveCodePro函数说明

| 函数名称                         | 函数描述                                                     | 参考文献 |
| -------------------------------- | ------------------------------------------------------------ | -------- |
| mAccCaliDescent                  | 基于加速度计矢量模相等的标定方法，代价函数值下降的搜索方法   | [1]      |
| mAccCaliHxJx1                    | 加速度计模型计算雅各比和海森矩阵                             |          |
| mAcce6PosCalibration             | 六位置修正加速度计标定法                                     | [2]      |
| GA算法框架（文件夹内所有函数）   | 遗传算法，可以用于非线性最优规划                             |          |
| ARkf                             | 马氏距离单因子抗差算法                                       | [4],[9]  |
| Adaptivekf156                    | 序贯量测+方差受限算法 15*6维，估计零偏                       |          |
| Adaptivekf186                    | 序贯量测+方差受限算法 18*6维 ，估计零偏，杆臂                |          |
| Adaptivekf196                    | 序贯量测+方差受限算法 18*6维 ，估计零偏，杆臂，同步时间      |          |
| BayesianFilter                   | 受污染的贝叶斯滤波算法                                       | [5]      |
| MRkf                             | 马氏距离单因子修正抗差算法                                   | [3]      |
| RobustKf                         | 多因子自适应抗差算法                                         |          |
| SH_KF156                         | Sage-Husa自适应卡尔曼滤波（序贯量测）                        |          |
| SlideWindow_Rkf                  | 基于滑动窗，马氏距离检测滤波算法                             |          |
| kf                               | 标准卡尔曼滤波                                               |          |
| kf0                              | 标准卡拉曼滤波，更准确的时间融合方式                         |          |
| kf_phi                           | 增加航向观测的kalman                                         |          |
| Virtual_Lever_Arm                | 虚拟杆臂法估计杆臂                                           | [6]      |
| kfCalibrated_Xu                  | 两步修正法标定陀螺                                           | [7]      |
| mGyro6PosCalibration             | 六位置法标定陀螺                                             |          |
| mfunRef43                        | 系统级标定法                                                 | [8]      |
| mgenRot                          | 生成旋转矩阵                                                 |          |
| NHC                              | 车辆不完整约束                                               |          |
| UKF_constrait                    | ukf用于车辆不完整约束                                        |          |
| Vehicle_Constraint               | 高程+航向+不完整约束                                         |          |
| AttAddKmtcCons                   | gnss outrage 解决方案: 姿态匹配+ 车体运动学速度匹配 20维传递对准算法 估计主子安装角，车体安装角，不能估计杆臂 |          |
| Att_Vel_AddKmtcCons              | gnss outrage 解决方案:姿态匹配+速度匹配  23维传递对准算法 估计主子安装角，杆臂，车体安装角 |          |
| Att_Vel_AddKmtcCons1             | gnss outrage 解决方案:增加正常段使用车辆约束                 |          |
| P_AdaptiveKf_NHC_PR              | 17维位置组合卡尔曼滤波 ，考虑GNSS缺失情况下使用运动学约束 , 模式识别（PR） |          |
| VmP_NHC                          | 车体速度辅助+位置，失锁时使用NHC                             |          |
| VmP_NHC_RTS                      | 车体速度辅助+位置，失锁时使用NHC/RTS                         | [11]     |
| VnP_NHC                          | 导航速度辅助+位置，失锁时使用NHC                             |          |
| ZeroVelocity_Kf                  | 零速修正kalman滤波                                           |          |
| ins_inpure                       | 纯惯性导航                                                   |          |
| RTS                              | RTS平滑核心代码                                              |          |
| smooth_RTS                       | RTS平滑，平滑区间与外部观测频率一直                          |          |
| smooth_RTS_all                   | RTS平滑，平滑区间取全段数据                                  |          |
| smooth_TKF                       | 双向平滑                                                     | [9]      |
| Aligni0                          | 粗对准                                                       |          |
| Arw_Vrw2std                      | ARW, VRW 与 噪声标准差的转换                                 |          |
| EarthParameter                   | 地球参数                                                     |          |
| Expand_axis_fill_figure          | 去除画图的空白边界                                           |          |
| KF_Phi                           | 组合导航离散线性模型                                         |          |
| Myavp2imu                        | 通过avp反向生成imu数据                                       |          |
| WaveDenoise                      | imu小波去噪                                                  |          |
| avp_update                       | 导航更新                                                     |          |
| feed_back_correct                | 反馈校正                                                     |          |
| fplot                            | 计算导航误差并画图                                           |          |
| genImuAsb                        | 产生非正交变换矩阵，从理想 b 系到 传感器系                   |          |
| globalParameter                  | 地球参数全局变量                                             |          |
| mToolLatLonErrorMeters           | 弧度转化为米                                                 |          |
| mxGetGravity                     | 计算重力值                                                   |          |
| reverseEarthParameter            | 反向滤波更新地球参数值                                       |          |
| rmse                             | 计算rmse                                                     |          |
| slideVarStd                      | 滑动计算数据的方差、标准差，均值                             |          |
| xkplot                           | 滤波估计值和协方差画图                                       |          |
| SINS153Vel_TransferAlignment     | 速度匹配15维传递对准算法 ，仅估计杆臂                        |          |
| SINS186Att_Vel_TransferAlignment | 姿态速度匹配18维传递对准算法 ，估计安装角和杆臂              |          |
| VPmaster                         | 主惯导速度位置匹配                                           |          |
| UKF156                           | 组合导航15维ukf                                              |          |
| UKFParameter                     | ukk参数设置                                                  |          |
| state_function                   | ukf的非线性模型                                              |          |
| ukf_filter                       | ukf滤波                                                      |          |
| utChange                         | ut变化                                                       |          |
| Wavelet_Transform（文件夹）      | 小波变换                                                     |          |
| SimuAcc                          | 加速度计仿真数据发生器                                       |          |
| TrjSim_INS                       | 轨迹发生器                                                   |          |

# NaveCodePro Demos

| 函数名称                  | 函数描述                     | 参考文献 |
| ------------------------- | ---------------------------- | -------- |
| Acc_Calibration_Sim_main  | 加速度计标定仿真，重复性测试 |          |
| Acc_Calibration_main      | 加速度计标定实测             |          |
| Gyro_Calibration_main     | 陀螺仪标定仿真，重复性测试   |          |
| Gyro_Calibration_sim_main | 陀螺仪标定实测               |          |
| Integrate_navi_real_main  | 抗差组合导航方法实测         |          |
| Integrate_navi_sim_main   | 抗差组合导航方法仿真         |          |
| Other_main                | NHC,UKF,平滑等函数测试       |          |
| Outrage_main              | GNSS失锁方案实际数据测试     |          |
| Transfer_Alignment_main   | 传递对准测试                 |          |
| TrjSim_main               | 产生仿真数据                 |          |
| data_process_main         | 数据处理                     |          |
| mainLeverArmEst           | 杆臂估计方案测试             |          |
| AllanCov（文件夹）        | Allan方差测试                |          |

# 参考文献

[1] Frosio I, Pedersini F, Borghese N A. Autocalibration of Triaxial MEMS Accelerometers With Automatic Sensor Model Selection[J]. IEEE Sensors Journal,2012,12(6):2100-2108.

[2] Xu T, Xu X, Xu D, et al. A Novel Calibration Method Using Six Positions for MEMS Triaxial Accelerometer[J]. IEEE Transactions on Instrumentation and Measurement, 2020, 70: 1-11. 

[3] Chen J , Shu-Bi Z . A Novel Adaptively-Robust Strategy Based on the Mahalanobis Distance for GPS/INS Integrated Navigation Systems[J]. Sensors, 2018, 18(3):695.

[4]. Chang G . Robust Kalman filtering based on Mahalanobis distance as outlier judging criterion[J]. Journal of Geodesy, 2014, 88(4):391-401.

[5]. 孙增析，邓志东, 一种对成片连续野值不敏感的鲁棒滤波. 清华大学学报（自然科学版）, 1994.

[6] Borko, A., I. Klein and G. Even-Tzur, GNSS/INS Fusion with Virtual Lever-Arm Measurements. Sensors, 2018. 18(7): p. 2228.

[7]邹泽兰，徐同旭，徐祥，赵鹤鸣. 基于两步修正法的MEMS三轴陀螺仪标定方法[j],仪器仪表学报

[8]Zhou Q, Yu G, Li H, et al. A Novel MEMS Gyroscope In-Self Calibration Approach[J]. Sensors, 2020, 20(18): 5430. 

[9]Liu H, Nassar S, El-Sheimy N. Two-filter smoothing for accurate ins/gps land-vehicle navigation in urban centers [J]. IEEE Transactions on Vehicular Technology, 2010, 59(9): 4256-4267

[10] Yang Y X, Cui X Q. Adaptively robust filter with multi adaptive factors [J]. Survey Review, 2008, 40(309): 260-270.

[11] Hang G ,  Guo J ,  Min Y , et al. A weighted combination filter with nonholonomic constrains for integrated navigation systems[J]. Advances in Space Research, 2015, 55(5):1470-1476.

# 致谢

NaveCodePro是在严恭敏老师的工具箱衍生而来，代码均与组合导航相关，code思路为文献复现和本人的原创设计。感谢课题组徐祥老师，徐同旭师兄对NavecodePro中部分代码的贡献。

# 版权说明

NaveCodePro不涉及也不允许任何商用，代码为本人对三年研究课题的整理，如有代码疑问或侵权行为，可联系邮箱：zelanzou@163.com.
