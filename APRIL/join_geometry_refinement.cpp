#include "join_geometry_refinement.h"

// topology masks according to de-9im
char insideCode[] = "T*F**F***";
// char trueInsideCode[] = "T*F*FF***";
char coveredbyCode1[] = "T*F**F***";
char coveredbyCode2[] = "*TF**F***";
char coveredbyCode3[] = "**FT*F***";
char coveredbyCode4[] = "**F*TF***";
char containsCode[] = "T*****FF*";
char coversCode1[] = "T*****FF*";
char coversCode2[] = "*T****FF*";
char coversCode3[] = "***T**FF*";
char coversCode4[] = "****T*FF*";
char meetCode1[] = "FT*******"; 
char meetCode2[] = "F**T*****"; 
char meetCode3[] = "F***T****"; 
char equalCode[] = "T*F**FFF*"; 
char disjointCode[] = "FF*FF****";
char intersectCode1[] = "T********";
char intersectCode2[] = "*T*******";
char intersectCode3[] = "***T*****";
char intersectCode4[] = "****T****";

//DE-9IM masks
//define topological masks for refinement

// a within b
boost::geometry::de9im::mask withinMask("T*F*FF***");
// a covered by b 
boost::geometry::de9im::mask coveredbyMask("**F*TF***");
// a and b meet
boost::geometry::de9im::mask meetMask1(meetCode1);
boost::geometry::de9im::mask meetMask2(meetCode2);
boost::geometry::de9im::mask meetMask3(meetCode3);
// a and b are equal
boost::geometry::de9im::mask equalMask(equalCode);
// a and b are disjoint
boost::geometry::de9im::mask disjointMask(disjointCode);
// a contains b
boost::geometry::de9im::mask containsMask(containsCode);



Polygon loadPolygonGeometry(uint &recID, unordered_map<uint,unsigned long> &offsetMap, ifstream &fin){
	Polygon pol(recID);
	int readID;
	int vertexCount, polygonCount;
	double x,y;

	double polxMin = numeric_limits<int>::max();
	double polyMin = numeric_limits<int>::max();
	double polxMax = -numeric_limits<int>::max();
	double polyMax = -numeric_limits<int>::max();

	//search the map for the specific polygon offset
	unordered_map<uint,unsigned long>::const_iterator got = offsetMap.find(recID);
	if(got != offsetMap.end()){
		//set read offset
		fin.seekg(got->second-fin.tellg(), fin.cur);		
		//read rec ID
		fin.read((char*) &readID, sizeof(int));
		//read vertex count
		fin.read((char*) &vertexCount, sizeof(int));

		pol.vertices.reserve(vertexCount);

		for(int i=0; i<vertexCount; i++){
			fin.read((char*) &x, sizeof(double));
			fin.read((char*) &y, sizeof(double));

			pol.vertices.emplace_back(x,y);

			//save polygon's min/max (for mbr)
			polxMin = min(polxMin, x);
			polyMin = min(polyMin, y);
			polxMax = max(polxMax, x);
			polyMax = max(polyMax, y);
		}

		//add MBR to polygon
		pol.mbr.pMin.x = polxMin;
		pol.mbr.pMin.y = polyMin;
		pol.mbr.pMax.x = polxMax;
		pol.mbr.pMax.y = polyMax;
	}

	return pol;
}

polygon loadPolygonGeometryBOOST(uint &recID, unordered_map<uint,unsigned long> &offsetMap, ifstream &fin){
	polygon pol;
	int readID;
	int vertexCount, polygonCount;
	double x,y;

	//search the map for the specific polygon offset
	unordered_map<uint,unsigned long>::const_iterator got = offsetMap.find(recID);
	if(got != offsetMap.end()){ 
		//set read offset
		fin.seekg(got->second-fin.tellg(), fin.cur);		
		//read rec ID
		fin.read((char*) &readID, sizeof(int));
		//read vertex count
		fin.read((char*) &vertexCount, sizeof(int));
		for(int i=0; i<vertexCount; i++){
			fin.read((char*) &x, sizeof(double));
			fin.read((char*) &y, sizeof(double));

			pol.outer().push_back(point_xy(x,y));
		}
	}

	boost::geometry::correct(pol);
	return pol;
}

linestring loadLinestringGeometryBOOST(uint &recID, unordered_map<uint,unsigned long> &offsetMap, ifstream &fin){
	linestring ls;
	int readID;
	int vertexCount;
	double x,y;

	//search the map for the specific polygon offset
	unordered_map<uint,unsigned long>::const_iterator got = offsetMap.find(recID);
	if(got != offsetMap.end()){ 
		//set read offset
		fin.seekg(got->second-fin.tellg(), fin.cur);
		//read rec ID
		fin.read((char*) &readID, sizeof(int));
		//read vertex count
		fin.read((char*) &vertexCount, sizeof(int));
		for(int i=0; i<vertexCount; i++){
			fin.read((char*) &x, sizeof(double));
			fin.read((char*) &y, sizeof(double));
			ls.emplace_back(x,y);
		}
	}
	return ls;
}


MBR getCMBR(Polygon &polA, Polygon &polB){
	MBR cmbr;
	cmbr.pMin.x = max(polA.mbr.pMin.x, polB.mbr.pMin.x);
	cmbr.pMin.y = max(polA.mbr.pMin.y, polB.mbr.pMin.y);
	cmbr.pMax.x = min(polA.mbr.pMax.x, polB.mbr.pMax.x);
	cmbr.pMax.y = min(polA.mbr.pMax.y, polB.mbr.pMax.y);
	return cmbr;
}

bool refinementWithIDsLinestring(uint &idA, uint &idB, unordered_map<uint,unsigned long> &offsetMapR, unordered_map<uint,unsigned long> &offsetMapS, ifstream &finR, ifstream &finS){
	/*  BOOST GEOMETRY REFINEMENT */
	polygon boostPolygonR = loadPolygonGeometryBOOST(idA, offsetMapR, finR);
	linestring boostLinestringS = loadLinestringGeometryBOOST(idB, offsetMapS, finS);
	refinementCandidatesR += boost::geometry::num_points(boostPolygonR);
	refinementCandidatesS += boost::geometry::num_points(boostLinestringS);

    return boost::geometry::intersects(boostPolygonR, boostLinestringS);
   
}

bool refinementWithIDs(uint &idA, uint &idB, unordered_map<uint,unsigned long> &offsetMapR, unordered_map<uint,unsigned long> &offsetMapS, ifstream &finR, ifstream &finS){
	/*  BOOST GEOMETRY REFINEMENT */
	polygon boostPolygonR = loadPolygonGeometryBOOST(idA, offsetMapR, finR);
	polygon boostPolygonS = loadPolygonGeometryBOOST(idB, offsetMapS, finS);
	refinementCandidatesR += boost::geometry::num_points(boostPolygonR);
	refinementCandidatesS += boost::geometry::num_points(boostPolygonS);

    return boost::geometry::intersects(boostPolygonR, boostPolygonS);
   
}

bool refinementWithinWithIDs(uint &idA, uint &idB, unordered_map<uint,unsigned long> &offsetMapR, unordered_map<uint,unsigned long> &offsetMapS, ifstream &finR, ifstream &finS){
	/*  BOOST GEOMETRY REFINEMENT */
	polygon boostPolygonR = loadPolygonGeometryBOOST(idA, offsetMapR, finR);
	polygon boostPolygonS = loadPolygonGeometryBOOST(idB, offsetMapS, finS);
	refinementCandidatesR += boost::geometry::num_points(boostPolygonR);
	refinementCandidatesS += boost::geometry::num_points(boostPolygonS);

    return boost::geometry::covered_by(boostPolygonR, boostPolygonS);
   
}

static inline std::string createMask(polygon PolygonR, polygon PolygonS) {
	boost::geometry::de9im::matrix matrix = boost::geometry::relation(PolygonR , PolygonS);
	return matrix.str();
}

static inline bool compareDe9imChars(char character, char char_mask) {
    if (character != 'F' && char_mask == 'T') {
        // character is 0,1,2 and char_mask is T
        return true;
    } else if (character == 'F' && char_mask == 'F'){
        // both are F
        return true;
    } else {
        // no match
        return false;
    }
}

static inline bool compareMasks(std::string &de9imCode, char* maskCode) {
    for(int i=0; i<9; i++) {
        if (de9imCode[i] == '*' || maskCode[i] == '*' || compareDe9imChars(de9imCode[i], maskCode[i])){
            continue;
        } else {
            return false;
        }
    }
    return true;
}

static int refineFindRelation(polygon polygonR, polygon polygonS){
        // get the mask codes
        std::string code = createMask(polygonR, polygonS);

        //disjoint
        if(compareMasks(code, disjointCode)){
            return DISJOINT;
        }

		// equals
        if(compareMasks(code, equalCode)){
            return EQUAL;
        }

		// R contains S
		if(compareMasks(code, containsCode)) {
			return R_CONTAINS_S;
		}

		// R contained by S
		if(compareMasks(code, insideCode)) {
			return S_CONTAINS_R; // R contained by S
		}
		
		// meet
        if(compareMasks(code, meetCode1) ||
                    compareMasks(code, meetCode2) ||
                    compareMasks(code, meetCode3)){
            return MEET;
        }

		return OVERLAP;

		



        /**
         * it definitely intersects at this point
        */
        /*
         // check equality first because it is a subset of covers and covered by
         if(compareMasks(code, equalCode)){
             return TR_EQUAL;
         }
         */

        // covers
        // if(compareMasks(code, coversCode1) || compareMasks(code, coversCode2) || compareMasks(code, coversCode3)|| compareMasks(code, coversCode4)){
        //     // first check contains because it is a subset of covers
        //     if(compareMasks(code, containsCode)){
        //         return TR_CONTAINS;
        //     }
        //     return TR_COVERS;
        // }

        
       /*
        // meet
        if(compareMasks(code, meetCode1) ||
                    compareMasks(code, meetCode2) ||
                    compareMasks(code, meetCode3)){
            return TR_MEET;
        }

        // else return intersects
        return TR_INTERSECT;
       */
        //return REFINEMENT;
    }

int refinement_DE9IM_WithIDs(uint &idA, uint &idB, unordered_map<uint,unsigned long> &offsetMapR, unordered_map<uint,unsigned long> &offsetMapS, ifstream &finR, ifstream &finS){
    /*  BOOST GEOMETRY REFINEMENT TO DETECT TOPOLOGICAL RELATIONSHIP */
    polygon boostPolygonR = loadPolygonGeometryBOOST(idA, offsetMapR, finR);
    polygon boostPolygonS = loadPolygonGeometryBOOST(idB, offsetMapS, finS);

    // //disjoint
    // if(boost::geometry::relate(boostPolygonR, boostPolygonS, disjointMask)){
    //     return DISJOINT;
    // }
 
    // //equal
    // if(boost::geometry::relate(boostPolygonR, boostPolygonS, equalMask)){
    //     return EQUAL;
    // }

	// // contains
	// if(boost::geometry::relate(boostPolygonR, boostPolygonS, containsMask)){
    //     return R_CONTAINS_S;
    // }
	// if(boost::geometry::relate(boostPolygonS, boostPolygonR, containsMask)){
    //     return S_CONTAINS_R;
    // }

    // //meet
    // if(boost::geometry::relate(boostPolygonR, boostPolygonS, meetMask1) ||
    //     boost::geometry::relate(boostPolygonR, boostPolygonS, meetMask2) ||
    //     boost::geometry::relate(boostPolygonR, boostPolygonS, meetMask3)){
    //     return MEET;
    // }
 
    // //else return overlap
    // return OVERLAP;

	return refineFindRelation(boostPolygonR, boostPolygonS);
 
}

