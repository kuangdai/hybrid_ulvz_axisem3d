//
//  InparamYAML.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 8/31/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  YAML parser for input parameters
//  based on mini-yaml:
//  https://github.com/jimmiebergmann/mini-yaml

#ifndef InparamYAML_hpp
#define InparamYAML_hpp

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Weverything"
#include "yaml/Yaml.hpp"
#pragma clang diagnostic pop

#include "bstring.hpp"
#include <map>

class InparamYAML {
public:
    // constructor
    InparamYAML(const std::string &name):
    mName(name), mRoot() {
        // nothing
    }
    
    // parse
    void parse(const std::string &fname);
    
    // verbose
    std::string verbose() const;
    
    // get value by keyword vector
    template <typename T>
    T getv(const std::vector<std::string> &keyword) const {
        // initialize nodePtr with root
        const Yaml::Node *nodePtr = &mRoot;
        
        // recursively increase depth of nodePtr
        std::string errKey = "";
        for (int depth = 0; depth < keyword.size(); depth++) {
            // update key series for error messege
            errKey += keyword[depth];
            
            // array or object?
            if (keyword[depth].front() == '[') {
                ////////////////////// array //////////////////////
                // check and get index
                if (keyword[depth].back() != ']') {
                    throw std::runtime_error("InparamYAML::getv || "
                                             "Error parsing array index. || "
                                             "Keyword: " + errKey + " || "
                                             "In YAML: " + mName);
                }
                const std::string &indexStr =
                keyword[depth].substr(1, keyword[depth].size() - 2);
                int index;
                try {
                    index = bstring::cast<int>(indexStr, "InparamYAML::getv");
                } catch(...) {
                    throw std::runtime_error("InparamYAML::getv || "
                                             "Error parsing array index. || "
                                             "Keyword: " + errKey + " || "
                                             "In YAML: " + mName);
                }
                
                // check and get array element
                if (nodePtr->IsSequence()) {
                    if (index < nodePtr->Size()) {
                        // access array element by iter because
                        // mini-yaml does not provide const random access
                        auto iter = nodePtr->Begin();
                        for (int imove = 0; imove < index; imove++) {
                            iter++;
                        }
                        // update valuePtr to array element
                        nodePtr = &((*iter).second);
                    } else {
                        throw std::runtime_error("InparamYAML::getv || "
                                                 "Index out of range. || "
                                                 "Keyword: " + errKey + " || "
                                                 "In YAML: " + mName);
                    }
                } else {
                    throw std::runtime_error("InparamYAML::getv || "
                                             "Non-array parent node. || "
                                             "Keyword: " + errKey + " || "
                                             "In YAML: " + mName);
                }
            } else {
                ////////////////////// object //////////////////////
                // check and get object member
                if (nodePtr->IsMap()) {
                    // access object member by iter because
                    // mini-yaml does not provide const random access
                    auto it = nodePtr->Begin();
                    for (; it != nodePtr->End(); it++) {
                        if ((*it).first == keyword[depth]) {
                            break;
                        }
                    }
                    if (it == nodePtr->End()) {
                        throw std::runtime_error("InparamYAML::getv || "
                                                 "Error finding keyword. || "
                                                 "Keyword: " + errKey + " || "
                                                 "In YAML: " + mName);
                    }
                    // update valuePtr to array element
                    nodePtr = &((*it).second);
                } else {
                    throw std::runtime_error("InparamYAML::getv || "
                                             "Non-object parent node. || "
                                             "Keyword: " + errKey + " || "
                                             "In YAML: " + mName);
                }
            }
            // update key series for error messege
            if (depth < keyword.size() - 1) {
                errKey += ":";
            }
        }
        
        // get string and use bstring to cast
        return bstring::cast<T>(nodePtr->As<std::string>(),
                                "InparamYAML::getv for keyword " + errKey);
    }
    
    // get value by keyword string
    template <typename T>
    T gets(const std::string &keyword) const {
        return getv<T>(bstring::split(keyword, ":"));
    }
    
    // get value with string-typed options
    template <typename T>
    T getsWithOptions(const std::string &keyword,
                      const std::map<std::string, T> &options) const {
        const std::string &res = gets<std::string>(keyword);
        try {
            return options.at(res);
        } catch (...) {
            return bstring::cast<T>(res, "InparamYAML::gets "
                                    "for keyword " + keyword);
        }
    }
    
    // get value with string-typed limits
    template <typename T>
    T getsWithLimits(const std::string &keyword,
                     const std::map<std::string, T> &limits) const {
        const std::string &res = gets<std::string>(keyword);
        try {
            return limits.at(res);
        } catch (...) {
            throw std::runtime_error("InparamYAML::gets || "
                                     "Unacceptable value: " + res + " || "
                                     "Keyword: " + keyword + " || "
                                     "In YAML: " + mName);
        }
    }
    
private:
    const std::string mName;
    Yaml::Node mRoot;
};

#endif /* InparamYAML_hpp */
