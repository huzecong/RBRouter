#include "HYPlotPS.h"
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <algorithm>
#include <cmath>
#include <map>
#include <cassert>
#include <limits>

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

        int r = (int)(arc.color_.r_*225);
        int g = (int)(arc.color_.g_*225);
        int b = (int)(arc.color_.b_*225);
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

void plot_postscript(string filename, const vector<Point> &point,
                     const RBRoutingPlan &plan) {
    double size = (plan.height + plan.width) / 2 / 25;

    vector<PSLine> lines;
    vector<PSArc> arcs;
    PSColor black;
    black.r_ = black.g_ = black.b_ = 0.0;
    black.alpha_ = 1.0;
    for (const Point &p : point) {
        PSArc arc;
        arc.centerX_ = p.x;
        arc.centerY_ = p.y;
        arc.radius_ = size / 4;
        arc.startAngle_ = 0.0;
        arc.endAngle_ = 2 * M_PI;
        arc.width_ = size / 4;
        arc.color_ = black;
        arcs.push_back(arc);
    }
    srand(time(0));
    for (int x = 0; x < plan.path.size(); ++x) {
        const vector<ID> &vec = plan.path[x];
        PSColor color;
        color.r_ = (float)(rand() % (12 + 1)) / 12;
        color.g_ = (float)(rand() % (12 + 1)) / 12;
        color.b_ = (float)(rand() % (12 + 1)) / 12;
        color.alpha_ = 1.0;
//		debug(color.r_ << " " << color.g_ << " " << color.b_);
        for (int i = 0; i + 1 < vec.size(); ++i) {
            Point a = point[vec[i]];
            Point b = point[vec[i + 1]];
            PSLine line;
            line.x1_ = a.x + ((float)x / plan.path.size() * 2 - 1) * size / 3;
            line.y1_ = a.y + ((float)x / plan.path.size() * 2 - 1) * size / 3;
            line.x2_ = b.x + ((float)x / plan.path.size() * 2 - 1) * size / 3;
            line.y2_ = b.y + ((float)x / plan.path.size() * 2 - 1) * size / 3;
            line.width_ = size / 3;
            line.color_ = color;
            lines.push_back(line);
        }
    }

    PSPlot plot;
    plot.setLines(lines);
    plot.setArcs(arcs);
    plot.draw(filename);
}