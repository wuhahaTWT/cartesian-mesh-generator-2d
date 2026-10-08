#pragma once
#include <filesystem>
#include <fstream>
#include <iterator>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>
#include <utility>

namespace cartmesh2d::cli {
// The native v1 case is deliberately a flat map of the CLI's typed options.
// It uses the same validation and defaults as a direct invocation. The desktop
// request archive is a different format and is not a path-bound native case.
class FlowCaseJson {
    std::string text_;
    std::size_t pos_=0;
    [[noreturn]] void fail(const std::string& message) const {
        throw std::invalid_argument("case-file JSON at byte "+std::to_string(pos_)+": "+message);
    }
    void space(){while(pos_<text_.size() && (text_[pos_]==' '||text_[pos_]=='\n'||text_[pos_]=='\r'||text_[pos_]=='\t'))++pos_;}
    bool take(char c){space();if(pos_<text_.size() && text_[pos_]==c){++pos_;return true;}return false;}
    void expect(char c){if(!take(c))fail(std::string("expected '")+c+"'");}
    unsigned hex(){
        unsigned value=0;
        for(int i=0;i<4;++i){
            if(pos_==text_.size())fail("incomplete Unicode escape");
            const char c=text_[pos_++];value*=16;
            if(c>='0'&&c<='9')value+=static_cast<unsigned>(c-'0');
            else if(c>='a'&&c<='f')value+=static_cast<unsigned>(c-'a'+10);
            else if(c>='A'&&c<='F')value+=static_cast<unsigned>(c-'A'+10);
            else fail("invalid Unicode escape");
        }
        return value;
    }
    std::string string(){
        expect('"');std::string result;
        while(pos_<text_.size()){
            const auto c=static_cast<unsigned char>(text_[pos_++]);
            if(c=='"')return result;
            if(c<32)fail("unescaped control character");
            if(c!='\\'){result+=static_cast<char>(c);continue;}
            if(pos_==text_.size())fail("incomplete escape");
            switch(text_[pos_++]){
            case '"':result+='"';break;case '\\':result+='\\';break;case '/':result+='/';break;
            case 'b':result+='\b';break;case 'f':result+='\f';break;
            case 'n':result+='\n';break;case 'r':result+='\r';break;case 't':result+='\t';break;
            case 'u':{
                auto cp=hex();
                if(cp>=0xd800&&cp<=0xdbff){
                    if(pos_+2>text_.size()||text_.substr(pos_,2)!="\\u")fail("missing low surrogate");
                    pos_+=2;const auto low=hex();
                    if(low<0xdc00||low>0xdfff)fail("invalid low surrogate");
                    cp=0x10000+(cp-0xd800)*1024+low-0xdc00;
                }else if(cp>=0xdc00&&cp<=0xdfff)fail("unpaired low surrogate");
                if(cp==0)fail("NUL is not a CLI value");
                if(cp<0x80)result+=static_cast<char>(cp);
                else if(cp<0x800){result+=static_cast<char>(0xc0|(cp>>6));result+=static_cast<char>(0x80|(cp&63));}
                else if(cp<0x10000){result+=static_cast<char>(0xe0|(cp>>12));result+=static_cast<char>(0x80|((cp>>6)&63));result+=static_cast<char>(0x80|(cp&63));}
                else{result+=static_cast<char>(0xf0|(cp>>18));result+=static_cast<char>(0x80|((cp>>12)&63));result+=static_cast<char>(0x80|((cp>>6)&63));result+=static_cast<char>(0x80|(cp&63));}
                break;
            }
            default:fail("invalid escape");
            }
        }
        fail("unterminated string");
    }
    std::string number(){
        space();const auto begin=pos_;
        if(pos_<text_.size()&&text_[pos_]=='-')++pos_;
        const auto digit=[&]{return pos_<text_.size()&&text_[pos_]>='0'&&text_[pos_]<='9';};
        if(!digit())fail("expected number");
        if(text_[pos_]=='0')++pos_;else while(digit())++pos_;
        if(pos_<text_.size()&&text_[pos_]=='.'){++pos_;if(!digit())fail("missing fractional digits");while(digit())++pos_;}
        if(pos_<text_.size()&&(text_[pos_]=='e'||text_[pos_]=='E')){
            ++pos_;if(pos_<text_.size()&&(text_[pos_]=='+'||text_[pos_]=='-'))++pos_;
            if(!digit())fail("missing exponent digits");
            while(digit())++pos_;
        }
        return text_.substr(begin,pos_-begin);
    }
public:
    static const std::set<std::string>& paths(){
        static const std::set<std::string> names{"mesh","output","boundary","export-boundaries","face-viscosity","initial-guess","initial-flux","restart"};return names;
    }
    explicit FlowCaseJson(std::string text):text_(std::move(text)){}
    std::map<std::string,std::string> read(){
        std::map<std::string,std::string> options;std::set<std::string> keys;std::string format;
        const std::set<std::string> strings{"case","linear-policy","convergence","coupling","steady-acceleration","viscous-stress","convection","outlet-backflow","pressure-preconditioner"};
        expect('{');
        if(!take('}'))do{
            const auto key=string();expect(':');if(!keys.insert(key).second)fail("duplicate key "+key);
            if(key=="format")format=string();
            else if(key=="options"){
                expect('{');if(!take('}'))do{
                    const auto name=string();expect(':');std::string value;
                    if(name.empty()||name.starts_with('-')||name=="case-file"||name=="help")fail("invalid option "+name);
                    if(name=="profile"){
                        space();if(text_.compare(pos_,4,"true")==0){value="true";pos_+=4;}
                        else if(text_.compare(pos_,5,"false")==0){value="false";pos_+=5;}
                        else fail("profile requires a boolean");
                    }else value=(strings.contains(name)||paths().contains(name))?string():number();
                    if(!options.emplace(name,value).second)fail("duplicate option "+name);
                    if(take('}'))break;
                    expect(',');
                }while(true);
            }else fail("unknown key "+key);
            if(take('}'))break;
            expect(',');
        }while(true);
        space();if(pos_!=text_.size())fail("trailing input");
        if(format!="cartmesh2d-native-flow-case-v1"||!keys.contains("options"))fail("expected native flow case v1 format and options");
        return options;
    }
};

inline std::vector<std::string> flowCaseArguments(int argc,char** argv){
    std::vector<std::string> explicitArgs;std::set<std::string> overrides;std::string file;
    for(int i=1;i<argc;++i){
        const std::string key=argv[i];
        if(key=="--case-file"){
            if(!file.empty()||i+1==argc)throw std::invalid_argument("provide exactly one --case-file path");
            file=argv[++i];continue;
        }
        explicitArgs.push_back(key);
        if(key.starts_with("--"))overrides.insert(key.substr(2));
        if(key!="--help"&&key!="--profile"&&i+1<argc)explicitArgs.emplace_back(argv[++i]);
    }
    if(file.empty())return explicitArgs;
    std::ifstream input(file);if(!input)throw std::runtime_error("cannot open case-file "+file);
    const std::string text{std::istreambuf_iterator<char>(input),{}};
    if(input.bad())throw std::runtime_error("cannot read case-file "+file);
    auto options=FlowCaseJson(text).read();std::vector<std::string> result;
    for(auto& [name,value]:options){
        if(overrides.contains(name))continue;
        if(name=="profile"){if(value=="true")result.push_back("--profile");continue;}
        if(FlowCaseJson::paths().contains(name)){
            if(value.empty())throw std::invalid_argument("empty case-file path: "+name);
            value=(std::filesystem::absolute(file).parent_path()/value).lexically_normal().string();
        }
        result.push_back("--"+name);result.push_back(value);
    }
    result.insert(result.end(),explicitArgs.begin(),explicitArgs.end());return result;
}
} // namespace cartmesh2d::cli
