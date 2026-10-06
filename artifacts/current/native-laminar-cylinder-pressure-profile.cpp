#include "cartmesh2d/io/MeshIO2D.hpp"
#include "cartmesh2d/fv/FvMesh2D.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

// Post-process the native wall-face pressure written by cartmesh2d_flow_cli.
// This is a profile/convergence diagnostic on accepted product fields, not an
// independent equation solve.  A constant pressure gauge is removed before
// profile comparisons; closed-body pressure force is gauge invariant.
namespace {

struct Sample {
    double angle{};
    double pressure{};
    double length{};
    double x{};
    double y{};
    double sx{};
    double sy{};
};

struct Profile {
    std::string label;
    std::vector<Sample> wall;
    std::vector<double> sampled;
    double perimeter{};
    double centreX{};
    double centreY{};
    double mean{};
    double minimum{std::numeric_limits<double>::infinity()};
    double maximum{-std::numeric_limits<double>::infinity()};
    double minimumAngle{};
    double maximumAngle{};
    double forceX{};
    double forceY{};
    double rms{};
    double totalVariation{};
};

std::vector<std::string> splitCsv(std::string line) {
    std::vector<std::string> result;
    std::size_t begin=0;
    while(true) {
        const auto comma=line.find(',',begin);
        result.push_back(line.substr(begin,comma==std::string::npos?comma:comma-begin));
        if(comma==std::string::npos)break;
        begin=comma+1;
    }
    return result;
}

double wrap(double angle) {
    const double twoPi=2*std::acos(-1.);
    angle=std::fmod(angle,twoPi);
    return angle<0?angle+twoPi:angle;
}

double interpolate(const std::vector<Sample>& samples,double angle) {
    angle=wrap(angle);
    const auto it=std::lower_bound(samples.begin(),samples.end(),angle,
        [](const Sample& sample,double value){return sample.angle<value;});
    const auto right=it==samples.end()?samples.begin():it;
    const auto left=right==samples.begin()?std::prev(samples.end()):std::prev(right);
    double a=left->angle,b=right->angle,x=angle;
    const double twoPi=2*std::acos(-1.);
    if(right==samples.begin())b+=twoPi;
    if(x<a)x+=twoPi;
    const double weight=(x-a)/(b-a);
    return (1-weight)*left->pressure+weight*right->pressure;
}

Profile readProfile(const std::string& label,const std::string& meshPath,
                    const std::string& facesPath,std::size_t sampleCount) {
    const auto read=cartmesh2d::readCm2dTopology(meshPath);
    if(!read.valid())throw std::runtime_error(read.error);
    const auto mesh=cartmesh2d::fv::makeFvMesh2D(read.topology);
    std::ifstream input(facesPath);
    if(!input)throw std::runtime_error("cannot read face field: "+facesPath);
    std::string line;
    if(!std::getline(input,line))throw std::runtime_error("empty face field");
    const auto header=splitCsv(line);
    auto column=[&](const std::string& name) {
        const auto it=std::find(header.begin(),header.end(),name);
        if(it==header.end())throw std::runtime_error("missing face column: "+name);
        return static_cast<std::size_t>(it-header.begin());
    };
    const auto faceColumn=column("face"),pressureColumn=column("pressure"),wallColumn=column("wall");
    Profile result;result.label=label;
    std::vector<bool> seen(mesh.faces.size());
    while(std::getline(input,line)) {
        const auto fields=splitCsv(line);
        const auto face=static_cast<std::size_t>(std::stoull(fields.at(faceColumn)));
        if(face>=mesh.faces.size()||seen[face])throw std::runtime_error("invalid duplicate face row");
        seen[face]=true;
        if(std::stoi(fields.at(wallColumn))==0)continue;
        const auto& geometry=mesh.faces[face];
        if(geometry.neighbour)throw std::runtime_error("wall marker on internal face");
        const double pressure=std::stod(fields.at(pressureColumn));
        const double length=std::hypot(geometry.areaVector.x,geometry.areaVector.y);
        result.wall.push_back({0.,pressure,length,geometry.centre.x,geometry.centre.y,
            geometry.areaVector.x,geometry.areaVector.y});
        result.perimeter+=length;
        result.centreX+=length*geometry.centre.x;
        result.centreY+=length*geometry.centre.y;
        result.mean+=length*pressure;
        result.forceX+=pressure*geometry.areaVector.x;
        result.forceY+=pressure*geometry.areaVector.y;
    }
    if(std::find(seen.begin(),seen.end(),false)!=seen.end())throw std::runtime_error("incomplete face field");
    if(result.wall.size()<8||!(result.perimeter>0))throw std::runtime_error("insufficient wall profile");
    // The benchmark boundary is a regular polygon.  Its perimeter centroid is
    // the circle centre even when product faces split each polygon segment.
    // Using the global origin would spuriously phase-shift a translated body.
    result.centreX/=result.perimeter;
    result.centreY/=result.perimeter;
    result.mean/=result.perimeter;
    for(auto& sample:result.wall)
        sample.angle=wrap(std::atan2(sample.y-result.centreY,sample.x-result.centreX));
    std::sort(result.wall.begin(),result.wall.end(),[](const auto& a,const auto& b){return a.angle<b.angle;});
    for(auto& sample:result.wall)sample.pressure-=result.mean;
    result.sampled.reserve(sampleCount);
    const double twoPi=2*std::acos(-1.);
    double square=0;
    for(std::size_t i=0;i<sampleCount;++i) {
        const double value=interpolate(result.wall,twoPi*i/sampleCount);
        result.sampled.push_back(value);
        if(value<result.minimum) {
            result.minimum=value;
            result.minimumAngle=twoPi*i/sampleCount;
        }
        if(value>result.maximum) {
            result.maximum=value;
            result.maximumAngle=twoPi*i/sampleCount;
        }
        square+=value*value;
    }
    result.rms=std::sqrt(square/sampleCount);
    for(std::size_t i=0;i<sampleCount;++i)
        result.totalVariation+=std::abs(result.sampled[(i+1)%sampleCount]-result.sampled[i]);
    return result;
}

void number(std::ostream& out,double value) {
    if(!std::isfinite(value))throw std::runtime_error("nonfinite diagnostic");
    out<<value;
}

} // namespace

int main(int argc,char** argv) {
    if(argc<7||(argc-1)%3)throw std::runtime_error(
        "usage: cylinder-pressure-profile label mesh.cm2d faces.csv [label mesh.cm2d faces.csv ...]");
    constexpr std::size_t sampleCount=16384;
    std::vector<Profile> profiles;
    for(int i=1;i<argc;i+=3)profiles.push_back(readProfile(argv[i],argv[i+1],argv[i+2],sampleCount));
    std::cout<<std::setprecision(17)<<"{\"angular_samples\":"<<sampleCount<<",\"cases\":[";
    for(std::size_t i=0;i<profiles.size();++i) {
        if(i)std::cout<<',';
        const auto& p=profiles[i];
        std::cout<<"{\"label\":\""<<p.label<<"\",\"wall_faces\":"<<p.wall.size()
                 <<",\"perimeter_m\":";number(std::cout,p.perimeter);
        std::cout<<",\"perimeter_centre_m\": [";number(std::cout,p.centreX);
        std::cout<<',';number(std::cout,p.centreY);std::cout<<']';
        std::cout<<",\"pressure_gauge_mean_m2_s2\":";number(std::cout,p.mean);
        std::cout<<",\"gauge_removed_profile_m2_s2\":{\"minimum\":";number(std::cout,p.minimum);
        std::cout<<",\"minimum_angle_rad\":";number(std::cout,p.minimumAngle);
        std::cout<<",\"maximum\":";number(std::cout,p.maximum);
        std::cout<<",\"maximum_angle_rad\":";number(std::cout,p.maximumAngle);
        std::cout<<",\"rms\":";number(std::cout,p.rms);
        std::cout<<",\"total_variation\":";number(std::cout,p.totalVariation);
        std::cout<<"},\"pressure_force_per_density_m4_s2\":{\"x\":";number(std::cout,p.forceX);
        std::cout<<",\"y\":";number(std::cout,p.forceY);std::cout<<"}}";
    }
    std::cout<<"],\"successive_profile_differences_m2_s2\":[";
    for(std::size_t i=1;i<profiles.size();++i) {
        if(i>1)std::cout<<',';
        const auto& a=profiles[i-1];const auto& b=profiles[i];
        double square=0,maximum=0,maximumAngle=0;
        for(std::size_t j=0;j<sampleCount;++j) {
            const double difference=b.sampled[j]-a.sampled[j];
            square+=difference*difference;
            if(std::abs(difference)>maximum) {
                maximum=std::abs(difference);
                maximumAngle=2*std::acos(-1.)*j/sampleCount;
            }
        }
        std::cout<<"{\"coarse\":\""<<a.label<<"\",\"fine\":\""<<b.label<<"\",\"rms\":";
        number(std::cout,std::sqrt(square/sampleCount));std::cout<<",\"maximum\":";
        number(std::cout,maximum);std::cout<<",\"maximum_angle_rad\":";
        number(std::cout,maximumAngle);std::cout<<'}';
    }
    std::cout<<"]}\n";
}
