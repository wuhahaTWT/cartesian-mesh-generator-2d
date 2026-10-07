#include "cartmesh2d/fv/detail/CompatibleFlowCheckpoint2D.hpp"
#include <bit>
#include <charconv>
#include <cmath>
#include <iomanip>
#include <istream>
#include <limits>
#include <ostream>
#include <sstream>
#include <string_view>
namespace cartmesh2d::fv {
namespace {
using Access=detail::CompatibleCheckpointAccess2D;
using Data=detail::CompatibleCheckpointData2D;
constexpr std::uint64_t offset=14695981039346656037ULL,prime=1099511628211ULL;
void hashWord(std::uint64_t& hash,std::uint64_t word) {
    for(int shift=56;shift>=0;shift-=8){hash^=(word>>shift)&255U;hash*=prime;}
}
[[noreturn]] void fail(const char* message){throw std::invalid_argument(std::string("Compatible checkpoint: ")+message);}
void token(std::istream& in,std::string_view expected) {
    std::string text;if(!(in>>text)||text!=expected)fail("unexpected token or truncated file");
}
std::uint64_t integer(std::istream& in,int base=10) {
    std::string text;if(!(in>>text))fail("missing integer");
    std::uint64_t value=0;const auto result=std::from_chars(text.data(),text.data()+text.size(),value,base);
    if(result.ec!=std::errc{}||result.ptr!=text.data()+text.size())fail("invalid integer");
    return value;
}
void real(double value){if(!std::isfinite(value))fail("nonfinite state or metadata");}
std::array<double,5> metrics(const CompatibleFlowMetrics2D& m) {
    return {m.cellMomentum,m.faceMomentum,m.divergence,m.stateChange,m.residualNorm};
}
void validate(const Data& data) {
    if(!data.context||!data.cells||!data.faces||!data.acceptedIterations||
       data.normalizedState.size()!=9*data.cells+4*data.faces)fail("invalid state extent or acceptance count");
    real(data.nextPseudoStep);if(data.nextPseudoStep<=0)fail("invalid next pseudo-step");
    for(double value:data.normalizedState)real(value);
    for(double value:metrics(data.metrics)){real(value);if(value<0)fail("negative residual metadata");}
}
void hexWord(std::ostream& out,std::uint64_t value) {
    std::array<char,16> buffer{};const auto result=std::to_chars(buffer.data(),buffer.data()+buffer.size(),value,16);
    const auto length=static_cast<std::size_t>(result.ptr-buffer.data());
    out<<std::string(16-length,'0');out.write(buffer.data(),static_cast<std::streamsize>(length));
}
}
std::size_t CompatibleFlowCheckpoint2D::acceptedIterations()const{return Access::get(*this).acceptedIterations;}
double CompatibleFlowCheckpoint2D::nextPseudoStep()const{return Access::get(*this).nextPseudoStep;}
CompatibleFlowMetrics2D CompatibleFlowCheckpoint2D::lastAcceptedMetrics()const{return Access::get(*this).metrics;}
void writeCompatibleFlowCheckpoint2D(std::ostream& out,const CompatibleFlowCheckpoint2D& checkpoint) {
    const auto& data=Access::get(checkpoint);validate(data);
    out<<"CARTMESH2D_COMPATIBLE_CHECKPOINT 1\nCOUNTS ";
    std::uint64_t hash=offset;
    for(auto n:{data.cells,data.faces,data.context->size(),data.normalizedState.size()}){
        out<<std::to_string(n)<<' ';hashWord(hash,n);
    }
    out<<"\nACCEPTED "<<std::to_string(data.acceptedIterations)<<'\n';hashWord(hash,data.acceptedIterations);
    const auto word=[&](std::uint64_t value){hexWord(out,value);out<<'\n';hashWord(hash,value);};
    out<<"NEXT_STEP\n";word(std::bit_cast<std::uint64_t>(data.nextPseudoStep));
    out<<"METRICS\n";for(double x:metrics(data.metrics))word(std::bit_cast<std::uint64_t>(x));
    out<<"CONTEXT\n";for(auto x:*data.context)word(x);
    out<<"STATE\n";for(double x:data.normalizedState)word(std::bit_cast<std::uint64_t>(x));
    out<<"CHECKSUM ";hexWord(out,hash);out<<"\nEND\n";
    if(!out)throw std::runtime_error("Compatible checkpoint: write failed");
}
CompatibleFlowCheckpoint2D readCompatibleFlowCheckpoint2D(std::istream& in,const FvMesh2D& mesh) {
    validateFvMesh2D(mesh);
    token(in,"CARTMESH2D_COMPATIBLE_CHECKPOINT");token(in,"1");token(in,"COUNTS");
    const std::array<std::size_t,4> expected{mesh.cells.size(),mesh.faces.size(),
        detail::compatibleCheckpointContextWords2D(mesh),9*mesh.cells.size()+4*mesh.faces.size()};
    std::uint64_t hash=offset;
    for(auto n:expected){const auto value=integer(in);if(value!=n)fail("mesh or state extent differs");hashWord(hash,value);}
    Data data;data.cells=expected[0];data.faces=expected[1];
    token(in,"ACCEPTED");const auto accepted=integer(in);
    if(accepted==0||accepted>std::numeric_limits<std::size_t>::max())fail("invalid acceptance count");
    data.acceptedIterations=static_cast<std::size_t>(accepted);hashWord(hash,accepted);
    const auto word=[&](){const auto value=integer(in,16);hashWord(hash,value);return value;};
    const auto number=[&](){const auto value=std::bit_cast<double>(word());real(value);return value;};
    token(in,"NEXT_STEP");data.nextPseudoStep=number();
    token(in,"METRICS");data.metrics={number(),number(),number(),number(),number()};
    token(in,"CONTEXT");auto context=std::make_shared<std::vector<std::uint64_t>>();context->reserve(expected[2]);
    for(std::size_t i=0;i<expected[2];++i)context->push_back(word());data.context=std::move(context);
    token(in,"STATE");data.normalizedState.reserve(expected[3]);
    for(std::size_t i=0;i<expected[3];++i)data.normalizedState.push_back(number());
    token(in,"CHECKSUM");if(integer(in,16)!=hash)fail("checksum mismatch");
    token(in,"END");std::string trailing;if(in>>trailing)fail("trailing data");if(in.bad())fail("read failed");
    validate(data);
    // Actual physical/discrete configuration is checked against freshly
    // prepared native equations by resumeCompatibleIncompressible2D.
    return Access::make(std::move(data));
}
}
