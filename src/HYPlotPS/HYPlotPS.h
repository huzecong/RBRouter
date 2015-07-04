#ifndef HY_PLOTPS_H
#define HY_PLOTPS_H

#include <vector>
#include <string>
#include <fstream>

struct PSColor {
    float r_;
    float g_;
    float b_;
    float alpha_;
};

struct PSArc {
    float centerX_;
    float centerY_;
    float radius_;
    float startAngle_;
    float endAngle_;
    float width_;
    PSColor color_;
};

struct PSLine {
    float x1_;
    float y1_;
    float x2_;
    float y2_;
    float width_;
    PSColor color_;
};

class PSPlot {
    std::vector<PSLine> lines_;
    std::vector<PSArc> arcs_;

    void plotCoords(std::ostream &out, float startX, float startY, float step, float blX, float blY, float rtX, float rtY);
    void plotGrids(std::ostream &out, float startX, float startY, float step, float blX, float blY, float rtX, float rtY);
    void plotArcs(std::ostream &out, float startX, float startY, float step, float blX, float blY);
    void plotLines(std::ostream &out, float startX, float startY, float step, float blX, float blY);
    void computeBBox(float &blX, float &blY, float &rtX, float &rtY);

    float computeX(float x, float startX, float blX, float step)
    {
        return (startX+(x-blX)*step);
    }

    float computeY(float y, float startY, float blY, float step)
    {
        return (startY+(y-blY)*step);
    }

public:
    void setLines(std::vector<PSLine> &lines) { lines_ = lines; }
    void setArcs(std::vector<PSArc> &arcs) { arcs_ = arcs; }
    void draw(std::string fileName);
};

#endif

