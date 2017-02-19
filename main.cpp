#include"CameraCalibration.h"
#include"StereoMatch.h"
#include<ctime>
Mat pointClouds;
static void onMouse(int event, int x, int y, int /*flags*/, void* /*param*/);
enum { BM, SGBM, VAR };
int main() {
	int APPROACH;
	StereoCalib sc;
	StereoMatch sm;
	StereoCalib::CornerDatas cornerDatas[2];
	StereoCalib::CameraParams cameraParams[2];
	StereoCalib::StereoParams stereoParams;
	StereoCalib::RemapMatrixs remapMatrixs, remapMatrixs_left, remapMatrixs_right;
	VideoCapture cam, cam_left, cam_right;
	Mat frame, empty, im_left, im_right, disparity, match_left, match_right, im_disparity,depth;
	empty.empty();
	//两个相机分别找棋盘角点
	/*for (int cameraID = 1;cameraID >=1;cameraID--) {
		cam.open(cameraID);
		bool first = true;
		int imageCount = 0;
		while (1) {
			cam>> frame;
			if (first) {
				sc.initCornerData(14, frame.size(), Size(6, 9), 36, cornerDatas[cameraID]);
				first = false;
			}
			int success = sc.detectCorners(frame, empty, cornerDatas[cameraID], imageCount);
			imshow("Findchessboardcorners", frame);
			if (success) {
				imageCount++;
				if (imageCount == 14)break;
				waitKey(0);
			}
			else if (waitKey(30) == 27)
				break;
		}
		cam.release();
	}
	destroyWindow("Findchessboardcorners");
	//sc.saveCornerData("cornerdata_left.yml", cornerDatas[0]);
	sc.saveCornerData("cornerdata_right.yml", cornerDatas[1]);
	//左相机定标
	//sc.calibrateSingleCamera(cornerDatas[0], cameraParams[0]);
	//sc.saveCameraParams(cameraParams[0], "cameraparams_left.yml");
	//右相机定标
	sc.calibrateSingleCamera(cornerDatas[1], cameraParams[1]);
	sc.saveCameraParams(cameraParams[1], "cameraparams_right.yml");
	//双目定标
	//cornerDatas[0].imagePoints2 = cornerDatas[1].imagePoints1;
	//sc.calibrateStereoCamera(cornerDatas[0], stereoParams,false);

	//sc.loadCornerData("cornerdata_left.yml", cornerDatas[0]);
	//sc.loadCornerData("cornerdata_right.yml", cornerDatas[1]);
	//sc.calibrateSingleCamera(cornerDatas[0], cameraParams[0]);
	//sc.calibrateSingleCamera(cornerDatas[1], cameraParams[1]);
	//cornerDatas[0].imagePoints2 = cornerDatas[1].imagePoints1;
	//sc.calibrateStereoCamera(cornerDatas[0], stereoParams, false);
	//sc.loadCameraParams("cameraparams_left.yml", cameraParams[0]);
	//sc.loadCameraParams("cameraparams_right.yml", cameraParams[1]);
	//单目矫正
	sc.rectifySingleCamera(cameraParams[0], remapMatrixs_left);
	sc.rectifySingleCamera(cameraParams[1], remapMatrixs_right);
	//双目矫正
	sc.rectifyStereoCamera(cornerDatas[0], stereoParams, remapMatrixs, StereoCalib::RECTIFY_BOUGUET);
	sc.saveCalibrationDatas("calibrationdatas.yml", StereoCalib::RECTIFY_BOUGUET, cornerDatas[0], stereoParams, remapMatrixs);*/

	//cam_left.open(2);
	//cam_right.open(0);
    cam_left.open("/Users/lrg/left.mp4");
    cam_right.open("/Users/lrg/right.mp4");
    //cam_left.set(CV_CAP_PROP_FRAME_WIDTH, 640);
    //cam_left.set(CV_CAP_PROP_FRAME_HEIGHT, 480);
    //cam_right.set(CV_CAP_PROP_FRAME_WIDTH, 640);
    //cam_right.set(CV_CAP_PROP_FRAME_HEIGHT, 480);
	cout << "重新定标？（y/n）:";
	char check;
	cin >> check;
	if (check == 'y') {
		bool first = true;
		int imageCount = 0;
		int chessborad_num, cols, rows, width;
		while (1) {
			cam_left >> im_left;
			cam_right >> im_right;
            if(im_left.empty()||im_right.empty())
                continue;
			if (first) {
				cout << "需检测棋盘数量：";
				cin >> chessborad_num;
				cout << "横向内角点数：";
				cin >> cols;
				cout << "纵向内角点数：";
				cin >> rows;
				cout << "棋盘格宽度（毫米）：";
				cin >> width;
				sc.initCornerData(chessborad_num/*20*/, im_left.size(), Size(cols, rows)/*6,9*/, width/*36*/, cornerDatas[0]);
				first = false;
			}
			int success = sc.detectCorners(im_left, im_right, cornerDatas[0], imageCount);
			imshow("Findchessboardcorners_left", im_left);
			imshow("Findchessboardcorners_right", im_right);
			if (success) {
				imageCount++;
				if (imageCount == chessborad_num)break;
				waitKey(0);
			}
			else if (waitKey(30) == 27)
				break;
		}
		sc.saveCornerData("/Users/lrg/Stereo/Stereo/cornerDatas.yml", cornerDatas[0]);
	}
	sc.loadCornerData("/Users/lrg/Stereo/Stereo/cornerDatas.yml", cornerDatas[0]);
	sc.calibrateStereoCamera(cornerDatas[0], stereoParams, true);
	sc.rectifyStereoCamera(cornerDatas[0], stereoParams, remapMatrixs, StereoCalib::RECTIFY_BOUGUET);
	sc.saveCalibrationDatas("/Users/lrg/Stereo/Stereo/calib_paras.xml", StereoCalib::RECTIFY_BOUGUET, cornerDatas[0], stereoParams, remapMatrixs);
	namedWindow("depth");
	setMouseCallback("depth", onMouse, 0);
	cout << "选择双目匹配算法：0-BM 1-SGBM 2-VAR :";
	cin >> APPROACH;
	while (1) {
		double minVal, maxVal;
		cam_left >> im_left;
		cam_right >> im_right;
        if(im_left.empty()||im_right.empty()){
            break;
        }
        resize(im_left, im_left, Size(640,480));
        resize(im_right, im_right, Size(640,480));
        clock_t start=clock();
        if(im_left.empty()||im_right.empty())
            continue;
		//sc.remapImage(im_left, im_right, im_left, im_right, remapMatrixs);
		//imshow("Rectified_left", im_left);
		//imshow("Rectified_right", im_right);
		//计算视差
		if (APPROACH == BM) {
			//int SADWindowSize = 9;
			//int numberOfDisparities = 16 * 8; /**< Range of disparity */
			//numberOfDisparities = numberOfDisparities > 0 ? numberOfDisparities : ((im_left.cols / 8) + 15) & -16;
			//sm.m_BM.state->roi1 = remapMatrixs.roi1;//左右视图的有效像素区域，一般由双目校正阶段的 cvStereoRectify 函数传递，也可以自行设定。
			//sm.m_BM.state->roi2 = remapMatrixs.roi2;//一旦在状态参数中设定了 roi1 和 roi2，OpenCV 会通过cvGetValidDisparityROI 函数计算出视差图的有效区域，在有效区域外的视差值将被清零。
			//sm.m_BM.state->preFilterSize=15;//预处理滤波器窗口大小,5-21,odd
			sm.m_BM.state->preFilterCap = 31; //63,1-31//预处理滤波器的截断值，预处理的输出值仅保留[-preFilterCap, preFilterCap]范围内的值,
			sm.m_BM.state->SADWindowSize = 19; //SAD窗口大小5-21
			sm.m_BM.state->minDisparity = 0; //64 最小视差，默认值为 0
			sm.m_BM.state->numberOfDisparities = ((640 / 8) + 15) & -16;//48; //128视差窗口，即最大视差值与最小视差值之差, 窗口大小必须是 16 的整数倍
			sm.m_BM.state->textureThreshold = 10;//低纹理区域的判断阈值。如果当前SAD窗口内所有邻居像素点的x导数绝对值之和小于指定阈值，则该窗口对应的像素点的视差值为 0
			sm.m_BM.state->uniquenessRatio = 15;//5-15 视差唯一性百分比， 视差窗口范围内最低代价是次低代价的(1 + uniquenessRatio/100)倍时，最低代价对应的视差值才是该像素点的视差，否则该像素点的视差为 0
			sm.m_BM.state->speckleWindowSize = 100;//检查视差连通区域变化度的窗口大小, 值为 0 时取消 speckle 检查
			sm.m_BM.state->speckleRange = 32;//视差变化阈值，当窗口内视差变化大于阈值时，该窗口内的视差清零
			sm.m_BM.state->disp12MaxDiff = 1;//左视差图（直接计算得出）和右视差图（通过cvValidateDisparity计算得出）之间的最大容许差异。超过该阈值的视差值将被清零。该参数默认为 -1，即不执行左右视差检查。
										//注意在程序调试阶段最好保持该值为 -1，以便查看不同视差窗口生成的视差效果。

										// 计算视差*/
			sm.bmMatch(im_left, im_right, disparity,match_left,match_right);
			//-- Check its extreme values
			//minMaxLoc(disparity, &minVal, &maxVal);
			//cout << "Min disp: Max value" << minVal << maxVal; //numberOfDisparities.= (maxVal - minVal)

															   //-- 4. Display it as a CV_8UC1 image
			//disparity.convertTo(disparity, CV_8U, 255 / (maxVal - minVal));//(numberOfDisparities*16.)
			//normalize(disparity, disparity, 0, 255, CV_MINMAX, CV_8UC1);    // obtain normalized image
			//GaussianBlur(disparity, disparity, Size(3, 3), 0);
			//imshow("disparity", disparity);
		}
		else if (APPROACH == SGBM) {
			sm.m_SGBM.preFilterCap = 63;
			int cn = 3;
			sm.m_SGBM.SADWindowSize = 9;
			sm.m_SGBM.P1 = 8 * cn*sm.m_SGBM.SADWindowSize*sm.m_SGBM.SADWindowSize;
			sm.m_SGBM.P2 = 32 * cn*sm.m_SGBM.SADWindowSize*sm.m_SGBM.SADWindowSize;
			sm.m_SGBM.minDisparity = 0;
			sm.m_SGBM.numberOfDisparities = ((640 / 8) + 15) & -16;//32;
			sm.m_SGBM.uniquenessRatio = 10;
			sm.m_SGBM.speckleWindowSize = 100;
			sm.m_SGBM.speckleRange = 32;
			sm.m_SGBM.disp12MaxDiff = 1;
			sm.m_SGBM.fullDP = 0;
			sm.sgbmMatch(im_left, im_right, disparity, match_left, match_right);
			//disparity.convertTo()
			//-- Check its extreme values
			//cv::minMaxLoc(disparity, &minVal, &maxVal);
			//cout << "Min disp: Max value" << minVal << maxVal; //numberOfDisparities.= (maxVal - minVal)

			//-- 4. Display it as a CV_8UC1 image
			//disparity.convertTo(disparity, CV_8U, 255 / (maxVal - minVal));//(numberOfDisparities*16.)
			//cv::normalize(disparity, disparity, 0, 255, CV_MINMAX, CV_8UC1);    // obtain normalized image
			//imshow("disparity", disparity);
		}
		else if (APPROACH == VAR) {
			int numberOfDisparities = 32;
			sm.m_VAR.levels = 3;                                 // ignored with USE_AUTO_PARAMS
			sm.m_VAR.pyrScale = 0.5;                             // ignored with USE_AUTO_PARAMS
			sm.m_VAR.nIt = 25;
			sm.m_VAR.minDisp = -(((640 / 8) + 15) & -16);
			sm.m_VAR.maxDisp = 0;
			sm.m_VAR.poly_n = 3;
			sm.m_VAR.poly_sigma = 0.0;
			sm.m_VAR.fi = 15.0f;
			sm.m_VAR.lambda = 0.03f;
			sm.m_VAR.penalization = sm.m_VAR.PENALIZATION_TICHONOV;   // ignored with USE_AUTO_PARAMS
			sm.m_VAR.cycle = sm.m_VAR.CYCLE_V;                        // ignored with USE_AUTO_PARAMS
			sm.m_VAR.flags = sm.m_VAR.USE_SMART_ID | sm.m_VAR.USE_AUTO_PARAMS | sm.m_VAR.USE_INITIAL_DISPARITY | sm.m_VAR.USE_MEDIAN_FILTERING;
			sm.varMatch(im_left, im_right, disparity, match_left, match_right);
		}
		//GaussianBlur(disparity, disparity, Size(11, 11), 0);
		sm.getPointClouds(disparity, pointClouds);
		//Mat topDownview;
		//sm.getTopDownView(pointClouds, topDownview);
		//sm.savePointClouds(pointClouds, "pointclouds.txt");
		sm.getDisparityImage(disparity, im_disparity,false);
		//获取深度
        double max_z = 1.8e4;
        for (int y = 0; y < pointClouds.rows; y++)
        {
            for (int x = 0; x < pointClouds.cols; x++)
            {
                cv::Vec3f point = pointClouds.at<cv::Vec3f>(y, x);
                if (fabs(point[2] - max_z) < FLT_EPSILON || point[2] > max_z || point[2]<0)
                    point[2]=0;
            }
        }
        vector<Mat> pos;
        split(pointClouds,pos);
        pos[2].convertTo(depth,CV_8U,255.0/1.8e4);
		imshow("Match_left", match_left);
		imshow("Match_right", match_right);
		imshow("disparity", im_disparity);
        imshow("depth",depth);
		//imshow("topdownview", topDownview);
		//imshow("pointclouds", pointClouds);
        clock_t end=clock();
        cout<<(double)(end-start)/CLOCKS_PER_SEC<<endl;
        if (waitKey(1) == 27){
            sm.savePointClouds(pointClouds, "/Users/lrg/Stereo/Stereo/pointclouds.txt");
            waitKey();
        }
	}
	return 0;
}
static void onMouse(int event, int x, int y, int /*flags*/, void* /*param*/) {

	if (event == CV_EVENT_LBUTTONDOWN)
	{
		/*// 根据深度阈值进行二值化处理
		double maxVal = 0, minVal = 0;
		cv::Mat depth2 = cv::Mat::zeros(depth.rows, depth.cols, CV_8UC1);
		cv::minMaxLoc(depth, &minVal, &maxVal);
		//double thrVal = minVal * 1.5;
		//threshold(depth, depthThresh, thrVal, 255, CV_THRESH_BINARY_INV);
		//depthThresh.convertTo(depthThresh, CV_8UC1);
		depth2.convertTo(depth2, CV_8U, 255 / (maxVal - minVal));//(numberOfDisparities*16.)
		normalize(depth2, depth2, 0, 255, CV_MINMAX, CV_8UC1);    // obtain normalized image*/
		cv::Vec3f point = pointClouds.at<cv::Vec3f>(y, x);
		cout << "("<< x << "," << y << ") distance:" << point[2] << endl;
	}
}
