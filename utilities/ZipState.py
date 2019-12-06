import zipfile
import os
import json
import shutil

#-------------------------------------------------------------------------------
#
#  Definition of the ZipState class.
#
#-------------------------------------------------------------------------------
class ZipState(zipfile.ZipFile):
    
#-------------------------------------------------------------------------------
#
#  ZipState constructor. This opens or creates a new zip file and sets the
#  inital state flags.
#
#-------------------------------------------------------------------------------
    def __init__(self, file_name, state):
        super(ZipState, self).__init__(file_name, state)
        if 'flags.json' not in self:
            self.set_state(state='unchanged')
#-------------------------------------------------------------------------------
#
#  ZipState in operator for simple query if a file exists in the archive.
#
#-------------------------------------------------------------------------------
    def __contains__(self, file):
        return file in super(ZipState, self).namelist()

#-------------------------------------------------------------------------------
#
#  ZipState write files. Writes a list of files to the zip file. Files may
#  already exist in the zip file so files are first split into files to write
#  and files to replace.
#
#-------------------------------------------------------------------------------
    def write(self, files):
#  The rest of the code assumes you are working with an array of file. If files
#  is a single instance, wrap it in a list.
        if not isinstance(files, list):
            files = [files]
        
#  Split files in the files to write and files to replace.
        write_files = []
        replace_files = []
        merge_files = []
        for file in files:
            if file in self:
                if file.endswith('.zip'):
                    merge_files.append(file)
                replace_files.append(file)
            else:
                write_files.append(file)

#  Merge files.
        self.merge(merge_files)

#  Replace existing file.
        self.replace(replace_files)

#  Write new files.
        for file in write_files:
            super(ZipState, self).write(file)

#-------------------------------------------------------------------------------
#
#  ZipState replace files. If a file is to be replaced, the entire zip file must
#  be recreated since they are immutable. Create a temp zip file and copy all
#  the contents. Then replace the existing file with the new one. Since this
#  is coping the entire file, replacement should be performed in batches.
#
#-------------------------------------------------------------------------------
    def replace(self, files):
#  Create a tempory zip file. Copy existing files and write replacement files.
        with zipfile.ZipFile('temp.zip', 'w') as zip_ref:
            for item in super(ZipState, self).infolist():
                if item.filename in files:
                    zip_ref.write(item.filename)
                else:
                    zip_ref.writestr(item, super(ZipState, self).read(item.filename))

#  Close the currently open zip file. Replace it with the temp file then reopen
#  the zip archive.
        super(ZipState, self).close()

        os.remove(self.filename)
        os.rename('temp.zip', self.filename)

        super(ZipState, self).__init__(self.filename, 'a')

#-------------------------------------------------------------------------------
#
#  ZipState remove files. If a file is to be removed, the entire zip file must
#  be recreated since they are immutable. Create a temp zip file and copy all
#  the contents except the removed files. Since this is coping the entire file,
#  removal should be performed in batches.
#
#-------------------------------------------------------------------------------
    def remove(self, files):
#  The rest of the code assumes you are working with an array of file. If files
#  is a single instance, wrap it in a list.
        if isinstance(files, str):
            files = [files]

#  Check all the files for the first existing file. If at least one file exists
#  begin the removal process.
        for file in files:
            if file in self:
                
#  A file exists in the zip archive. Copy all the contents unless it is in the
#  replacement list.
                with zipfile.ZipFile('temp.zip', 'w') as zip_ref:
                    for item in super(ZipState, self).infolist():
                        if item.filename not in files:
                            zip_ref.writestr(item, super(ZipState, self).read(item.filename))

#  Close the currently open zip file. Replace it with the temp file then reopen
#  the zip archive.
                super(ZipState, self).close()

                os.remove(self.filename)
                os.rename('temp.zip', self.filename)
    
                super(ZipState, self).__init__(self.filename, 'a')
                break

#-------------------------------------------------------------------------------
#
#  ZipState merge files. If a file does not exist in the new file add to the new
#  file.
#
#-------------------------------------------------------------------------------
    def merge(self, files):
#  Make a temporary sub directory and change the working directory to that
#  directory.
        os.mkdir('temp_dir')
        os.chdir('temp_dir')

        for file in files:
#  Rename the new file to a temp and extract the old archive.
            super(ZipState, self).extract(file)

#  Extract all files from the new archive and write them to the old archive.
            with zipfile.ZipFile('../{}'.format(file), 'a') as new_zip_ref:
                new_zip_ref.extractall()
                with ZipState(file, 'a') as old_zip_ref:
                    old_zip_ref.write(new_zip_ref.namelist())
                
#  Clean up all the files extracted from the new file and remove the temp file.
                for file2 in new_zip_ref.namelist():
                    os.remove(file2)

#  Copy the current file back to the parent directory.
            shutil.copy2(file, '../')

#  Change working directory back to the parent.
        os.chdir('../')
        shutil.rmtree('temp_dir')

#-------------------------------------------------------------------------------
#
#  ZipState get the state dictionary. The state is store in a json file.
#
#-------------------------------------------------------------------------------
    def get_state(self):
        if 'flags.json' in super(ZipState, self).namelist():
            super(ZipState, self).extract('flags.json')
            with open('flags.json', 'r') as flags_file:
                return json.load(flags_file)
        else:
            return {}

#-------------------------------------------------------------------------------
#
#  ZipState set flags in the state dictionary. The state is store in a json
#  file.
#
#-------------------------------------------------------------------------------
    def set_state(self, **keywords):
        flags = self.get_state()
        
        for key, value in keywords.items():
            flags[key] = value

        with open('flags.json', 'w') as flags_file:
            json.dump(flags, flags_file)
        self.write('flags.json')
