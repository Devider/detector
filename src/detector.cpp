#include <opencv2/opencv.hpp>
#include <stdio.h>
#include <sys/time.h>
#include <iostream>



using namespace std;
using namespace cv;

struct CartesianLine {
	float k;
	float b;
};

struct Circle {
	Point center;
	int radius;
};

void points(Point* p1, Point* p2, double rho, double theta) {
	double a = cos(theta), b = sin(theta);
	double x0 = a * rho, y0 = b * rho;
	p1->x = cvRound(x0 + 1000 * (-b));
	p1->y = cvRound(y0 + 1000 * (a));
	p2->x = cvRound(x0 - 1000 * (-b));
	p2->y = cvRound(y0 - 1000 * (a));
	if (p1->x == p2->x)
		p1->x += 1;
	if (p1->y == p2->y)
		p1->y += 1;
}

void transfotmToCartesianLine(vector<CartesianLine>& toLines, vector<Vec2f>& fromLines) {
	toLines.clear();
	for (uint i = 0; i < fromLines.size(); i++){
		Vec2f& fromLine = fromLines[i];
		CartesianLine toLine;
		float rho = fromLine[0], theta = fromLine[1];
		Point p1, p2;
		points(&p1, &p2, rho, theta);
		if ((p2.x - p1.x) == 0)
			return;
		float b = ((p2.x * p1.y) - (p2.y * p1.x)) / (p2.x - p1.x);
		float k = (p1.y - b) / p1.x;
		toLine.b = b;
		toLine.k = k;
		toLines.push_back(toLine);
	}
}

void average (Point& p, vector<Point>& points){
	if (points.size() == 0)
		return;
	int x = 0,y = 0;
	for (uint i = 0; i < points.size();i++) {
		x += points[i].x;
		y += points[i].y;
	}
	x = x / points.size();
	y = y / points.size();
	p.x = x;
	p.y = y;
}

bool correct (int value, int min, int max){
	return value > min && value < max;
}

void addPoint(Point& p, vector<Point>& points){
	if (correct(p.x, 0, 640) && correct (p.y, 0, 480))
				points.push_back(p);
};

void detectIntersection(int currentPosition, vector<Point>& points, vector<CartesianLine>& lines) {
	int cp = currentPosition + 1;
	int k = 2;
	for (uint i = cp; i < lines.size(); i++){
		if (pow(lines[currentPosition].k - lines[i].k, 2.0) > k) {
			float x = (lines[i].b - lines[currentPosition].b) / (lines[currentPosition].k - lines[i].k);
			float y = lines[currentPosition].k * x + lines[currentPosition].b;
			Point p = Point(x, y);
			addPoint(p , points);
		}
	}
}

void detectIntersection(vector<Point>& points, vector<CartesianLine>& cLines) {
	if (cLines.size() < 2)
		return;
	for (uint i = 0; i < cLines.size(); i++){
		detectIntersection(i, points, cLines);
	}

}

bool detectIntersections(Mat& gray, vector<Point>& intersections) {
	bool result = false;
	vector<Vec2f> lines;
	HoughLines(gray, lines, 1, CV_PI / 180, 170, 0, 0);
	for (size_t i = 0; i < lines.size(); i++) {
		float rho = lines[i][0], theta = lines[i][1];
		Point p1, p2;
		points(&p1, &p2, rho, theta);
	}

	vector<CartesianLine> cLines(lines.size());

	transfotmToCartesianLine(cLines, lines);

	detectIntersection (intersections, cLines);
	Point intersection;
	average(intersection, intersections);
	return result;
}

void detectCircles(Mat& gray, vector<Circle>& circles) {
	vector<Vec3f> points;
	HoughCircles(gray, points, CV_HOUGH_GRADIENT, 2, gray.rows / 4, 200, 100, 200, 480);
	for (size_t i = 0; i < points.size(); i++) {
		Circle c;
		c.center = Point(cvRound(points[i][0]), cvRound(points[i][1]));
		c.radius = cvRound(points[i][2]);
		circles.push_back(c);
	}
}

void drawPoints(Mat& mat, vector<Point>& points, int radius, Scalar color){
	for (uint i = 0; i < points.size(); i++) {
		circle(mat, points[i], 3, color, -1, 8, 0);
	}
}

void drawPoints(Mat& mat, vector<Circle>& circles, int radius, Scalar color){
	vector<Point> points(circles.size());
	for (uint i = 0; i < circles.size(); i++) {
		points[i] = circles[i].center;
	}
	drawPoints(mat, points, radius, color);
}

int main(int argc, char** argv) {
	VideoCapture cap;
	cap.set(CV_CAP_PROP_FPS , 10);
	if (argc > 1)
		cap.open(string(argv[1]));
	else
		cap.open(0);

	if (!cap.isOpened())
		cout << "NOT OPENED!" << endl;

	Mat frame, edged;

	time_t start, stop;



	for (;;) {

		struct timeval  tv;
		gettimeofday(&tv, NULL);

		start = (tv.tv_sec) * 1000 + (tv.tv_usec) / 1000 ;

		vector<Point> intersections;
		vector<Circle> circles;
		cap >> frame;
//		Mat src = frame;

		if ( frame.empty() ){
			cout << "Frame is empty!" << endl;
			continue;
		}

		cvtColor( frame, frame, CV_RGB2GRAY );
		threshold( frame, frame, 97, 255, 0 );

		Canny(frame, edged, 50, 100);

		detectIntersections(edged, intersections);

		cout << "Lines intersections found: " << intersections.size() << "\n";

		Point p;
		average(p, intersections);



		detectCircles(edged, circles);

//		drawPoints(src, intersections, 3, Scalar(30, 50, 255));
//		circle(src, p, 5, Scalar(251,222,0), 5);
//		drawPoints(src, circles, 3, Scalar(250, 250, 0));

//		imshow("detected lines", edged);
//		imshow("circles", frame);
//		imshow("src", src);

		gettimeofday(&tv, NULL);
		stop = (tv.tv_sec) * 1000 + (tv.tv_usec) / 1000 ;

		cout << stop - start << "ms. passed\n";

		if (waitKey(30) >= 0)
			break;

	}
}

