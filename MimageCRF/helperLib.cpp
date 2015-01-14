#include "helperLib.h"
#include <vector>
#include <algorithm>

using std::vector;
using std::sort;


int readDirectories(std::string &dirNamestr, std::vector<string>& dirList, std::vector<string>& indirName, bool testDir ) {

      const char *dirName = dirNamestr.c_str ();

	  /* check command line arguments */
	  if (dirName == '\0') {
	    fprintf (stderr, "invalid directory");
	    return EXIT_FAILURE;
	  }

	  DIR *dir;
	  struct dirent *ent;
	  string name, fullPath;
	  std::vector<string> files;

	  /* open directory stream */
	  dir = opendir (dirName);
	  if (dir != NULL) {
		  /* print all the files and directories within directory */
	      while ((ent = readdir (dir)) != NULL) {
			  name = ent->d_name ;
			  if(name == "..")
				  continue;
			  else if(name == ".")
				  continue;
			  else if( std::find(name.begin(), name.end(), '.') != name.end()){
				  // this is a file because it is of the format X.Y
				  // have to make sure directories dont have a dot in their names
			  } else{ // assume the other has to be directory

				  // use this later // DirEntry->d_type == isFile
				  dirList.push_back(name);
			  }
	      }
		sort(dirList.begin(), dirList.end());

		if (testDir == 1) {

			int flag = 0;

			for(std::vector<string>::iterator b = dirList.begin(); b < dirList.end(); ++b){
				string foldn = *b;
				if (foldn.find("test") == 0) {
					flag = 1;
					break;
				}
			}

			if (flag == 0)
				throw CError("Directory does not contain directory: test* \n");

/*
			std::vector<string>::iterator it = std::find(dirList.begin(), dirList.end(), "test");
			if(it == dirList.end() ){
				throw CError("Directory does not contain directory: test \n");
			}
*/
		}

		for(std::vector<string>::iterator b = dirList.begin(); b < dirList.end(); ++b){
			  fullPath = dirName;
			  if( fullPath[fullPath.length() -1] != dirs)
				  fullPath += dirs;
			  fullPath += *b;
			  indirName.push_back(fullPath);

		}


		  closedir (dir);
	  } else {
	      /* could not open directory */
	      perror ("");
	      return EXIT_FAILURE;
	  }
	  return EXIT_SUCCESS;
}

// Read image from a directory.
int readImages(char *dirName, std::vector<CByteImage>& ims, bool testImage){
  /* check command line arguments */
  if (dirName == '\0') {
    fprintf (stderr, "invalid directory");
    return EXIT_FAILURE;
  }

  DIR *dir;
  struct dirent *ent;
  string name, fullPath;
  std::vector<string> files;

  /* open directory stream */
  dir = opendir (dirName);
  if (dir != NULL) {
	  /* print all the files and directories within directory */
      while ((ent = readdir (dir)) != NULL) {
		  name = ent->d_name ;
		  if(name == "..")
			  continue;
		  else if(name == ".")
			  continue;
		  else if( std::find(name.begin(), name.end(), '.') != name.end()){
			  fullPath = dirName;
			  if( fullPath[fullPath.length() -1] != dirs)
				  fullPath += dirs;
			  fullPath += name;
			  files.push_back(fullPath);
		  } else{ // a directory or unknow format.
		  }
      }
	sort(files.begin(), files.end());

	string testString = dirName;
	testString += dirs;
	testString += "test.png";
	if(testImage == 1){
		std::vector<string>::iterator it = std::find(files.begin(), files.end(), testString);
		if(it == files.end() ){
			throw CError("Directory does not contain test.png image \n");
		}
		CByteImage imTemp;
		ReadImage(imTemp, it->c_str());
		ims.push_back(imTemp);
		files.erase(it);
	}

	for(std::vector<string>::iterator b = files.begin(); b < files.end(); ++b){
		cout << *b << "\n";
		CByteImage imTemp;
		ReadImage(imTemp, b->c_str());
		ims.push_back(imTemp);
	}


	  closedir (dir);
  } else {
      /* could not open directory */
      perror ("");
      return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}


int readImages(char *dirName, std::vector<CImage2>& ims, bool testImage){
  /* check command line arguments */
  if (dirName == '\0') {
    fprintf (stderr, "invalid directory");
    return EXIT_FAILURE;
  }

  DIR *dir;
  struct dirent *ent;
  string name, fullPath;
  std::vector<string> files;

  /* open directory stream */
  dir = opendir (dirName);
  if (dir != NULL) {
	  /* print all the files and directories within directory */
      while ((ent = readdir (dir)) != NULL) {
		  name = ent->d_name ;
		  if(name == "..")
			  continue;
		  else if(name == ".")
			  continue;
		  else if( std::find(name.begin(), name.end(), '.') != name.end()){
			  fullPath = dirName;
			  if( fullPath[fullPath.length() -1] != dirs)
				  fullPath += dirs;
			  fullPath += name;
			  files.push_back(fullPath);
		  } else{ // a directory or unknow format.
		  }
      }
	sort(files.begin(), files.end());

	string testString = dirName;
	testString += dirs;
	testString += "test.png";
	if(testImage == 1){
		std::vector<string>::iterator it = std::find(files.begin(), files.end(), testString);
		if(it == files.end() ){
			throw CError("Directory does not contain test.png image \n");
		}
		CImage2 imTemp;
		char fname[100];
		strcpy(fname, it->c_str());
		read_png_file(imTemp, fname);
		ims.push_back(imTemp);
		files.erase(it);
	}

	for(std::vector<string>::iterator b = files.begin(); b < files.end(); ++b){
		cout << *b << "\n";
		CImage2 imTemp;
		char fname[100];
		strcpy(fname, b->c_str());
		read_png_file(imTemp, fname);
		ims.push_back(imTemp);
	}


	  closedir (dir);
  } else {
      /* could not open directory */
      perror ("");
      return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}


int readImages(string dirName, std::vector<CImage2>& ims){
  /* check command line arguments */
  if (dirName.length()==0) {
    fprintf (stderr, "invalid directory");
    return EXIT_FAILURE;
  }

  DIR *dir;
  struct dirent *ent;
  string name, fullPath;
  std::vector<string> files;

  /* open directory stream */
  dir = opendir (dirName.c_str());
  if (dir != NULL) {
	  /* print all the files and directories within directory */
      while ((ent = readdir (dir)) != NULL) {
		  name = ent->d_name ;
		  if(name == "..")
			  continue;
		  else if(name == ".")
			  continue;
		  else if (name == "Thumbs.db")
			  continue;
		  else if (name == "thumbs.db")
			  continue;
		  else if (name == ".DS_Store")
			  continue;
		  else if( std::find(name.begin(), name.end(), '.') != name.end()){
			  fullPath = dirName;
			  if( fullPath[fullPath.length() -1] != dirs)
				  fullPath += dirs;
			  fullPath += name;
			  files.push_back(fullPath);
		  } else{ // a directory or unknow format.
		  }
      }
	sort(files.begin(), files.end());

	for(std::vector<string>::iterator b = files.begin(); b < files.end(); ++b){
		cout << *b << "\n";
		CImage2 imTemp;
		char fname[100];
		strcpy(fname, b->c_str());
		read_png_file(imTemp, fname);
		ims.push_back(imTemp);
	}


	  closedir (dir);
  } else {
      /* could not open directory */
      perror ("");
      return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}




int readImages(string dirName, std::vector<CByteImage>& ims){
  /* check command line arguments */
  if (dirName.length()==0) {
    fprintf (stderr, "invalid directory");
    return EXIT_FAILURE;
  }

  DIR *dir;
  struct dirent *ent;
  string name, fullPath;
  std::vector<string> files;

  /* open directory stream */
  dir = opendir (dirName.c_str());
  if (dir != NULL) {
	  /* print all the files and directories within directory */
      while ((ent = readdir (dir)) != NULL) {
		  name = ent->d_name ;
		  if(name == "..")
			  continue;
		  else if(name == ".")
			  continue;
		  else if (name == "Thumbs.db")
			  continue;
		  else if (name == "thumbs.db")
			  continue;
		  else if (name == ".DS_Store")
			  continue;
		  else if( std::find(name.begin(), name.end(), '.') != name.end()){
			  fullPath = dirName;
			  if( fullPath[fullPath.length() -1] != dirs)
				  fullPath += dirs;
			  fullPath += name;
			  files.push_back(fullPath);
		  } else{ // a directory or unknow format.
		  }
      }
	sort(files.begin(), files.end());

	for(std::vector<string>::iterator b = files.begin(); b < files.end(); ++b){
		cout << *b << "\n";
		CByteImage imTemp;
		ReadImage(imTemp, b->c_str());
		ims.push_back(imTemp);
	}


	  closedir (dir);
  } else {
      /* could not open directory */
      perror ("");
      return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}



int readImages(string dirName, int subdir, std::vector<CImage2>& ims){
  /* check command line arguments */
  if (dirName.length()==0) {
    fprintf (stderr, "invalid directory");
    return EXIT_FAILURE;
  }

  DIR *dir;
  struct dirent *ent;
  string name, fullPath;
  std::vector<string> files;

  char dname[100];

  sprintf(dname, "%s/%d/", dirName.c_str(), subdir);

  string s2(dname);

  dirName = s2;

  /* open directory stream */
  dir = opendir (dname);
  if (dir != NULL) {
	  /* print all the files and directories within directory */
      while ((ent = readdir (dir)) != NULL) {
		  name = ent->d_name ;
		  if(name == "..")
			  continue;
		  else if(name == ".")
			  continue;
		  else if (name == "Thumbs.db")
			  continue;
		  else if (name == "thumbs.db")
			  continue;
		  else if (name == ".DS_Store")
			  continue;
		  else if( std::find(name.begin(), name.end(), '.') != name.end()){
			  fullPath = dirName;
			  if( fullPath[fullPath.length() -1] != dirs)
				  fullPath += dirs;
			  fullPath += name;
			  files.push_back(fullPath);
		  } else{ // a directory or unknow format.
		  }
      }
	sort(files.begin(), files.end());

	for(std::vector<string>::iterator b = files.begin(); b < files.end(); ++b){
		cout << *b << "\n";
		CImage2 imTemp;
		char fname[100];
		strcpy(fname, b->c_str());
		read_png_file(imTemp, fname);
		ims.push_back(imTemp);
	}


	  closedir (dir);
  } else {
      /* could not open directory */
      perror ("");
      return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}


int readImagesPaths(string dirName, std::vector<string>& filepaths){
  /* check command line arguments */
  if (dirName.length()==0) {
    fprintf (stderr, "invalid directory");
    return EXIT_FAILURE;
  }

  DIR *dir;
  struct dirent *ent;
  string name, fullPath;

  /* open directory stream */
  dir = opendir (dirName.c_str());
  if (dir != NULL) {
    /* print all the files and directories within directory */
      while ((ent = readdir (dir)) != NULL) {
      name = ent->d_name ;
      if(name == "..")
        continue;
      else if(name == ".")
        continue;
      else if (name == "Thumbs.db")
        continue;
      else if (name == "thumbs.db")
        continue;
      else if (name == ".DS_Store")
        continue;
      else if( std::find(name.begin(), name.end(), '.') != name.end()){
        fullPath = dirName;
        if( fullPath[fullPath.length() -1] != dirs)
          fullPath += dirs;
        fullPath += name;
        filepaths.push_back(fullPath);
      } else{ // a directory or unknow format.
      }
      }
  sort(filepaths.begin(), filepaths.end());

    closedir (dir);
  } else {
      /* could not open directory */
      perror ("");
      return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}




