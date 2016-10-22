function [nameWithoutExtension extensionOut bTest] = RemoveExtensionFromImageName(filename,bWarnings)
    if nargin<2
        bWarnings = 1;
    end
    extensions(1).name = '.vtk';
    extensions(2).name = '.VTK';
    extensions(3).name = '.nii';
    extensions(4).name = '.NII';
    extensions(5).name = '.gipl';
    extensions(6).name = '.GIPL';
    extensions(7).name = '.nii.gz';
    extensions(8).name = '.bm';
    extensions(9).name = '.raw';
    extensions(10).name = '.mha';
    extensions(11).name = '.MHA';
    extensions(12).name = '.cvi42wsx';
    extensions(13).name = '.par';
    extensions(14).name = '.rec';
    extensions(15).name = '.dcm';
    iTest  = 0;
    bTest = 0;
    nTests = numel(extensions);
    while (~bTest)&&(iTest<nTests)
        iTest = iTest+1;
        extension = extensions(iTest).name;
        [bTest,nameWithoutExtension,extensionOut] = TestExtension(filename,extension);
    end
    if(bWarnings)
        if bTest==0 && ~strcmp(filename,'none')
            fprintf('WARNING! In RemoveExtensionFromImageName: file type in image file ''%s'' not recognised! Extensions supported so far are:\n',filename);
            fprintf('  %s,',extensions(:).name);
            fprintf('\n');
        end
    end

function [bTest,nameWithoutExtension,extensionOut] = TestExtension(filename,extension)
    nameWithoutExtension = [];
    extensionOut = [];
    nChar = numel(extension);
    if numel(filename)<nChar
        bTest = 0;
        return;
    end
    if strcmp(filename(end-nChar+1:end),extension)
        nameWithoutExtension = filename(1:end-nChar);
        extensionOut = extension;
        bTest = 1;
    else
        bTest = 0;
    end