import zipfile
import os

class ZipState(zipfile.ZipFile):
    def __init__(self, file_name, state):
        super(ZipState, self).__init__(file_name, state)

    def write(self, files):
        replace = False
        
        if isinstance(files, basestring):
            files = [files]
        
        for file in super(ZipState, self).namelist():
            if file in files:
                replace = True

        if replace:
            self.replace(files)

        for file in files:
            super(ZipState, self).write(file)

    def replace(self, files):
        with zipfile.ZipFile('temp.zip', 'w') as zip_ref:
            for item in super(ZipState, self).infolist():
                if item.filename in files:
                    zip_ref.write(item.filename)
                    files.remove(item.filename)
                else:
                    zip_ref.writestr(item, super(ZipState, self).read(item.filename))

        super(ZipState, self).close()

        os.remove(self.filename)
        os.rename('temp.zip', self.filename)

        super(ZipState, self).__init__(self.filename, 'a')

    def remove(self, files):
        remove = False

        if isinstance(files, basestring):
            files = [files]

        for file in super(ZipState, self).namelist():
            if file in files:
                remove = True

        if remove:
            with zipfile.ZipFile('temp.zip', 'w') as zip_ref:
                for item in super(ZipState, self).infolist():
                    if item.filename not in files:
                        zip_ref.writestr(item, super(ZipState, self).read(item.filename))

            super(ZipState, self).close()

            os.remove(self.filename)
            os.rename('temp.zip', self.filename)
    
            super(ZipState, self).__init__(self.filename, 'a')
