#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include "HYPlotPS.h"

using namespace std;

void PSPlot::plotCoords(ostream &out, float startX, float startY, float step, float blX, float blY, float rtX, float rtY)
{
    out << "/Times-Roman findfont" << endl;

    out << std::max((int)step, 1) << " scalefont" << endl;
    out << "setfont" << endl;
    int x = 0, y = 0;
    for (;;)
    {
        if (x % 5 == 0)
        {
            out << std::max((float)(startX+step*x-0.2*step), 0.0f) << " " << std::max((float)(startY-step), 0.0f) << " moveto" << endl;
            out << "(" << blX+x << ") show" << endl; 
        }
        if ((++x) > (rtX-blX))
            break;
    }
    for (;;)
    {
        out << std::max((float)(startX-2*step), 0.0f) << " " << std::max((float)(startY+step*y-0.2*step), 0.0f) << " moveto" << endl;
        out << "(" << blY+y << ") show" << endl; 
        if ((++y) > (rtY-blY))
            break;
    }
}
 
void PSPlot::plotGrids(ostream &out, float startX, float startY, float step, float blX, float blY, float rtX, float rtY)
{
    out << "/verticals" << endl;
    out << "{ newpath" << endl;
    char line[1024];
    sprintf(line, "%.2f %.2f %.2f\n", startX, step, startX+(rtX-blX)*step);
    out << line;
    out << "{ " << startY << " moveto" << endl;
    sprintf(line, "%.2f %.2f rlineto } for\n", 0.0f, (rtY-blY)*step);
    out << line;
    out << "stroke" << endl;
    out << "} def" << endl;

    out << "/horizontals" << endl;
    out << "{ newpath" << endl;
    sprintf(line, "%.2f %.2f %.2f\n", startY, step, startY+(rtY-blY)*step);
    out << line;
    out << "{ " << startX << " exch moveto" << endl;
    sprintf(line, "%.2f %.2f rlineto } for\n", (rtX-blX)*step, 0.0f);
    out << line;
    out << "stroke" << endl;
    out << "} def" << endl;

    out << "[3 3] 0 setdash" << endl;
    out << "1 setlinejoin" << endl;
    sprintf(line, "%.2f setlinewidth\n", step/10.0);
    out << line;
    out << "verticals" << endl;
    out << "horizontals" << endl;
}

void PSPlot::plotArcs(ostream &out, float startX, float startY, float step, float blX, float blY)
{
    char line[1024];
    
    for (int i = 0; i < arcs_.size(); i++)
    {
        PSArc &arc = arcs_[i];

//        int r = (int)(arc.color_.r_*225);
//        int g = (int)(arc.color_.g_*225);
//        int b = (int)(arc.color_.b_*225);
		double r = arc.color_.r_;
		double g = arc.color_.g_;
		double b = arc.color_.b_;
        out << r << " " << g << " " << b << " setrgbcolor" << endl;
        float x = computeX(arcs_[i].centerX_, startX, blX, step);
        float y = computeY(arcs_[i].centerY_, startY, blY, step);
        sprintf(line, "%.2f setlinewidth\n", step/3.0*arc.width_);
        out << line;
        sprintf(line, "%.2f %.2f %.2f %.2f %.2f doArc\n", x, y, arc.radius_*step, arc.startAngle_*180.0/M_PI, arc.endAngle_*180.0/M_PI);
        out << line;
    }
}

void PSPlot::plotLines(ostream &out, float startX, float startY, float step, float blX, float blY)
{
    char lineStr[1024];
    out << "[] 0 setdash" << endl;
    out << "1 setlinejoin" << endl;

    //int r = 0x00;
    //int g = 0x99;
    //int b = 0x66;
    for (int i = 0; i < (int)lines_.size(); i++)
    {
        PSLine &line = lines_[i];
//        int r = (int)(line.color_.r_*225);
//        int g = (int)(line.color_.g_*225);
//        int b = (int)(line.color_.b_*225);
        double r = line.color_.r_;
        double g = line.color_.g_;
        double b = line.color_.b_;
        sprintf(lineStr, "%.2lf %.2lf %.2lf setrgbcolor\n", r, g, b);
        out << lineStr;
//        out << r << " " << g << " " << b << " setrgbcolor" << endl;
        float x1 = computeX(lines_[i].x1_, startX, blX, step);
        float x2 = computeX(lines_[i].x2_, startX, blX, step);
        float y1 = computeY(lines_[i].y1_, startY, blY, step);
        float y2 = computeY(lines_[i].y2_, startY, blY, step);
        sprintf(lineStr, "%.2f setlinewidth\n", step/3.0*line.width_);
        out << lineStr;
        out << "newpath" << endl;
        sprintf(lineStr, "%.2f %.2f %.2f %.2f doLine\n", x1, y1, x2, y2);
        out << lineStr;
    }
    out << "stroke" << endl;
}

void PSPlot::computeBBox(float &blX, float &blY, float &rtX, float &rtY)
{
    blX = numeric_limits<float>::max();
    blY = numeric_limits<float>::max();
    rtX = numeric_limits<float>::min();
    rtY = numeric_limits<float>::min();
    for (int i = 0; i < (int)lines_.size(); i++)
    {
        PSLine &line = lines_[i];
        if (blX > line.x1_)
            blX = line.x1_;
        if (blX > line.x2_)
            blX = line.x2_;
        if (blY > line.y1_)
            blY = line.y1_;
        if (blY > line.y2_)
            blY = line.y2_;
        if (rtX < line.x1_)
            rtX = line.x1_;
        if (rtX < line.x2_)
            rtX = line.x2_;
        if (rtY < line.y1_)
            rtY = line.y1_;
        if (rtY < line.y2_)
            rtY = line.y2_;
    }
    for (int i = 0; i < (int)arcs_.size(); i++)
    {
        PSArc &arc = arcs_[i];
        if (blX > arc.centerX_ - arc.radius_)
            blX = arc.centerX_ - arc.radius_;
        if (blY > arc.centerY_ - arc.radius_)
            blY = arc.centerY_ - arc.radius_;
        if (rtX < arc.centerX_ + arc.radius_)
            rtX = arc.centerX_ + arc.radius_;
        if (rtY < arc.centerY_ + arc.radius_)
            rtY = arc.centerY_ + arc.radius_;
    }
}

void PSPlot::draw(std::string fileName)
{
    ofstream out(fileName.c_str());

    float blX, blY, rtX, rtY;
    computeBBox(blX, blY, rtX, rtY);
    float pagewidth = 8.5 * 72;
    float pageheight = 11 * 72;
    float step = ((int)(std::min((pagewidth-100.0)/(rtX-blX), (pageheight-100.0)/(rtY-blY))*100))/100.0;
    float startX = ((int)((pagewidth-step*(rtX-blX))/2.0*100))/100.0;
    float startY = ((int)((pageheight-step*(rtY-blY))/2.0*100))/100.0;
    plotGrids(out, startX, startY, step, blX, blY, rtX, rtY);
    plotCoords(out, startX, startY, step, blX, blY, rtX, rtY);

    char line[1024];
    out << "/doCircle" << endl;
    out << "{ /radius exch def" << endl;
    out << "/ypos exch def" << endl;
    out << "/xpos exch def" << endl;
    out << "xpos ypos radius 0 360 arc fill stroke} def" << endl;

    out << "/doArc" << endl;
    out << "{ /endDeg exch def" << endl;
    out << "/startDeg exch def" << endl;
    out << "/radius exch def" << endl;
    out << "/ypos exch def" << endl;
    out << "/xpos exch def" << endl;
    out << "xpos ypos radius startDeg endDeg arc stroke} def" << endl;

    out << "/doBox" << endl;
    out << "{ /y2 exch def" << endl;
    out << "/x2 exch def" << endl;
    out << "/y1 exch def" << endl;
    out << "/x1 exch def" << endl;
    out << "x1 y1 moveto" << endl;
    out << "x1 y2 lineto" << endl;
    out << "x2 y2 lineto" << endl;
    out << "x2 y1 lineto" << endl;
    out << "closepath" << endl;
    out << "fill" << endl;
    out << "stroke} def" << endl;

    out << "/doLine" << endl;
    out << "{ /y2 exch def" << endl;
    out << "/x2 exch def" << endl;
    out << "/y1 exch def" << endl;
    out << "/x1 exch def" << endl;
    out << "x1 y1 moveto" << endl;
    out << "x2 y2 lineto" << endl;
    out << "stroke} def" << endl;

    plotLines(out, startX, startY, step, blX, blY);
    plotArcs(out, startX, startY, step, blX, blY);

    out << "showpage" << endl;
    out.close();
}

PSColor rand_color() {
    PSColor color;
    color.r_ = (float)(rand() % (12 + 1)) / 12;
    color.g_ = (float)(rand() % (12 + 1)) / 12;
    color.b_ = (float)(rand() % (12 + 1)) / 12;
    color.alpha_ = 1.0;
    return color;
}

struct PSCircle {
    Point p;
    double radius;
    PSCircle() {}
    PSCircle(Point _p, double _r) : p(_p), radius(_r) {}
    Point point_at_angle(double angle) const {
        return Point(p.x + radius * cos(angle), p.y + radius * sin(angle));
    }
};

void tangent(PSCircle A, PSCircle B, pair<double, double> *deg) {
	bool swapped = false;
    if (A.radius < B.radius) {
		swap(A, B);
		swapped = true;
	}
    double d = dist(A.p, B.p);
    double base = atan2(B.p.y - A.p.y, B.p.x - A.p.x);
    double angle = acos((A.radius - B.radius) / d);
    deg[0] = make_pair(base + angle, base + angle);
    deg[1] = make_pair(base - angle, base - angle);
    angle = acos((A.radius + B.radius) / d);
    deg[2] = make_pair(base + angle, M_PI + base + angle);
    deg[3] = make_pair(base - angle, M_PI + base - angle);
	if (swapped) {
		for (int i = 0; i < 4; ++i)
			swap(deg[i].first, deg[i].second);
	}
}

double epsilon;

void add_line(const Point &a, const Point &b, const PSColor &color, vector<PSLine> &lines) {
	PSLine line;
	line.x1_ = (float)a.x;
	line.y1_ = (float)a.y;
	line.x2_ = (float)b.x;
	line.y2_ = (float)b.y;
	line.width_ = (float)(epsilon / 2);
	line.color_ = color;
	lines.push_back(line);
}

void add_arc(const Point &p, double r, double s, double e,
			 const PSColor &color, vector<PSArc> &arcs) {
	if (s < 0) s += 2 * M_PI;
	if (e < 0) e += 2 * M_PI;
	if (s > e) {
		add_arc(p, r, min(s, 2 * M_PI), 2 * M_PI, color, arcs);
		add_arc(p, r, 0.0, e, color, arcs);
	} else {
		PSArc arc;
		arc.centerX_ = (float) p.x;
		arc.centerY_ = (float) p.y;
		arc.radius_ = (float) r;
		arc.startAngle_ = (float) s;
		arc.endAngle_ = (float) e;
		arc.width_ = (float) (epsilon / 2);
		arc.color_ = color;
		arcs.push_back(arc);
	}
};

// b is to the left of a
inline bool to_left_of(double a, double b) {
	a = rectify(a);
	b = rectify(b);
	if (b < a) b += 2 * M_PI;
	return a + M_PI > b - eps;
}

void plot_postscript(string filename, const vector<Point> &point,
                     const RBRoutingPlan &plan) {
	typedef RBRoutingPlan::Side Side;
    int max_layer = 1;
    for (const auto &x : plan.layer)
        for (auto y : x)
            max_layer = max(max_layer, (int)y);
    epsilon = min(plan.width, plan.height);
    for (int i = 0; i < point.size(); ++i)
        for (int j = i + 1; j < point.size(); ++j)
            epsilon = min(epsilon, dist(point[i], point[j]));
    epsilon /= max_layer * 5.0;
	epsilon = max(epsilon, min(plan.width, plan.height) / 500.0);

    vector<PSLine> lines;
    vector<PSArc> arcs;
	auto choose = [&](double x, Side dir, double *angle) {
		if (dir == RBRoutingPlan::LeftSide) {
			if (to_left_of(x, angle[0])) return angle[0];
			else return angle[1];
		} else {
			if (to_left_of(x, angle[0])) return angle[1];
			else return angle[0];
		}
	};
	auto draw = [&](const Point &a, const Point &b, double r,
					Side dir, const PSColor &color) {
		pair<double, double> _angle[4];
		PSCircle A(a, 0.0), B(b, r);
		tangent(A, B, _angle);	// only use first two
		double con = atan2(a.y - b.y, a.x - b.x);
		double angle[2];
		for (int i = 0; i < 2; ++i)
			angle[i] = _angle[i].second;
		double sel = choose(con, (RBRoutingPlan::Side) -dir, angle);
		Point p = B.point_at_angle(sel);
		add_line(a, p, color, lines);
		return sel;
	};

    PSColor black;
    black.r_ = black.g_ = black.b_ = 0.0;
    black.alpha_ = 1.0;
    for (const Point &p : point)
		add_arc(p, epsilon / 2, 0.0, 2 * M_PI, black, arcs);

    srand(time(0));
    set<PSColor> used_color;
    for (int x = 0; x < plan.path.size(); ++x) {
        const vector<ID> &vec = plan.path[x];
		print(plan.path[x]);
		print(plan.layer[x]);
        PSColor color;
        double luminance = 0.0;
        do {
            color = rand_color();
            luminance = (max(color.r_, max(color.g_, color.b_))
                         - min(color.r_, min(color.g_, color.b_))) / 2.0;
        } while (used_color.find(color) != used_color.end() || luminance > 8.0);
        used_color.insert(color);

		if (vec.size() == 2) {
			add_line(point[vec[0]], point[vec[1]], color, lines);
			continue;
		}
		double last_angle;
		last_angle = draw(point[vec[0]], point[vec[1]], plan.layer[x][1] *
														epsilon,
						  plan.direction[x][1], color);
        for (int i = 1; i + 2 < vec.size(); ++i) {
			pair<double, double> _angle[4], sel;
			const Point &a = point[vec[i]], &b = point[vec[i + 1]];
			PSCircle A(a, plan.layer[x][i] * epsilon);
			PSCircle B(b, plan.layer[x][i + 1] * epsilon);
			tangent(A, B, _angle);
			RBRoutingPlan::Side da = plan.direction[x][i];
			RBRoutingPlan::Side db = plan.direction[x][i + 1];
			if (da != db) {
				for (int i = 0; i < 2; ++i)
					swap(_angle[i], _angle[i + 2]);
			}
			double con = atan2(b.y - a.y, b.x - a.x);
			double angle[2];
			for (int i = 0; i < 2; ++i)
				angle[i] = _angle[i].first;
			double sa = choose(con, da, angle);
			for (int i = 0; i < 2; ++i)
				angle[i] = _angle[i].second;
			double sb = choose(con + M_PI, (Side) -db, angle);
			if (da == RBRoutingPlan::RightSide) {
				add_arc(a, plan.layer[x][i] * epsilon, last_angle, sa, color,
						arcs);
			} else {
				add_arc(a, plan.layer[x][i] * epsilon, sa, last_angle, color,
						arcs);
			}
			add_line(A.point_at_angle(sa), B.point_at_angle(sb), color, lines);
			last_angle = sb;
        }
		int len = vec.size() - 1;
		double cur_angle = draw(point[vec[len]], point[vec[len - 1]],
								plan.layer[x][len - 1] * epsilon,
								(Side) -plan.direction[x][len - 1], color);
		if (plan.direction[x][len - 1] == RBRoutingPlan::RightSide) {
			add_arc(point[vec[len - 1]], plan.layer[x][len - 1] * epsilon,
					last_angle, cur_angle, color, arcs);
		} else {
			add_arc(point[vec[len - 1]], plan.layer[x][len - 1] * epsilon,
					cur_angle, last_angle, color, arcs);
		}
    }

    PSPlot plot;
    plot.setLines(lines);
    plot.setArcs(arcs);
    plot.draw(filename);
}