function pmml = translate(self)
% Translate. Translate hyperparameters to valid GPR PMML
%   Return PMML as a string
    rootTag = 'PMML';
    xmlns = 'http://www.dmg.org/PMML-4_3';
    version='4.3';
    document = com.mathworks.xml.XMLUtils.createDocument(rootTag);
    pmmlElement = document.getDocumentElement;
    pmmlElement.setAttribute('xmlns',xmlns);
    pmmlElement.setAttribute('version',version);
    writeHeader(self,document,pmmlElement);
    writeDataDictionary(self,document,pmmlElement);
    writeModel(self,document,pmmlElement);
    pmml = xmlwrite(document);
end


function writeHeader(self,doc,parent)
    % Write the <header> tag
    header = doc.createElement('Header');
    header.setAttribute('copyright','DMG.org');
    parent.appendChild(header);
end


function writeDataDictionary(self,document,parent)
    % Write the <datadictionary> section
    nfields = size(self.xTrain,2) + size(self.yTrain,2);
    datadictionary = document.createElement('DataDictionary');
    datadictionary.setAttribute('numberOfFields',num2str(nfields));
    for i=1:size(self.xTrain,2)
        % Write the xTrain fields
        name = sprintf('x%i',i);
        datafield = document.createElement('DataField');
        datafield.setAttribute('dataType','double'); %TODO: Set this dynamically
        datafield.setAttribute('name',name); % TODO Custom field names
        datafield.setAttribute('optype','continuous'); % TODO: set dynamically
        datadictionary.appendChild(datafield);
    end
    for i=1:size(self.yTrain,2)
        % Write the yTrain fields
        name = sprintf('y%i',i);
        datafield = document.createElement('DataField');
        datafield.setAttribute('dataType','double'); %TODO: Set this dynamically
        datafield.setAttribute('name',name); % TODO: Custom field names
        datafield.setAttribute('optype','continuous'); %TODO Todo set dynamically
        datadictionary.appendChild(datafield);
    end
    parent.appendChild(datadictionary);
end


function writeModel(self,document,parent)
    % Write the <gaussianprocessmodel> section
    gpm = document.createElement('GaussianProcessModel');
    gpm.setAttribute('functionName','regression');
    gpm.setAttribute('modelName','Gaussian Process Model');
    writeMiningSchema(self,document,gpm);
    writeOutputSection(self,document,gpm);
    writeKernelSection(self,document,gpm);
    writeTrainingInstances(self,document,gpm);
    parent.appendChild(gpm);
end


function writeMiningSchema(self,document,parent)
    % Write the <miningschema> section
    miningschema = document.createElement('MiningSchema');
    for i=1:size(self.xTrain,2)
        defaultName = sprintf('x%i',i);
        miningfield = document.createElement('MiningField');
        miningfield.setAttribute('name',defaultName); % TODO Track actual field names
        miningfield.setAttribute('usageType','active');
        miningschema.appendChild(miningfield);
    end
    for i=1:size(self.yTrain,2)
        defaultName = sprintf('y%i',i);
        miningfield = document.createElement('MiningField');
        miningfield.setAttribute('name',defaultName); % TODO Track actual field names
        miningfield.setAttribute('usageType','predicted');
        miningschema.appendChild(miningfield);
    end
    parent.appendChild(miningschema);
end


function writeOutputSection(self,document,parent)
% Write the <output> tag
    outputfields = {'MeanValue','StandardDeviation'};
    output = document.createElement('Output');
    for i=1:length(outputfields)
        outputfield = document.createElement('OutputField');
        outputfield.setAttribute('dataType','double');
        outputfield.setAttribute('feature','predictedValue');
        outputfield.setAttribute('name',outputfields{i});
        outputfield.setAttribute('optype','continuous');
        output.appendChild(outputfield);
    end
    parent.appendChild(output);
end


function writeKernelSection(self,document,parent)
    % Write the <gaussianprocessmodel>
    kernel = num2str(self.covFunc);
    noise = exp(self.hyp.lik)^2;
    gamma = exp(self.hyp.cov(end))^2;
    lambda = exp(self.hyp.cov(1:end-1));
    lambdaArray = sprintf(' %f',lambda);
    lambdaArray = lambdaArray(2:end);

    kernelNode = document.createElement(kernel);
    lambdaNode = document.createElement('Lambda');
    arrayNode = document.createElement('array');

    kernelNode.setAttribute('gamma',num2str(gamma));
    kernelNode.setAttribute('noiseVariance',num2str(noise)); % Todao set NV term
    arrayNode.setAttribute('n',num2str(length(lambda)));
    arrayNode.setAttribute('type','real');
    arrayNode.setTextContent(lambdaArray);

    lambdaNode.appendChild(arrayNode);
    kernelNode.appendChild(lambdaNode);
    parent.appendChild(kernelNode);
end


function writeTrainingInstances(self,document,parent)
    % Write the <traininginstances> wrapper
    % This section wraps the instanceField (description) and inlinetable (data)
    nRows = size(self.xTrain,1);
    nCols = size(self.xTrain,2)+size(self.yTrain,2);
    trainingInstances = document.createElement('TrainingInstances');
    trainingInstances.setAttribute('recordCount',num2str(nRows));
    trainingInstances.setAttribute('fieldCount',num2str(nCols));
    trainingInstances.setAttribute('isTransformed','false');
    % Write nested fields
    writeInstanceFields(self,document,trainingInstances);
    writeInlineTable(self,document,trainingInstances);
    parent.appendChild(trainingInstances);
end


function writeInstanceFields(self,document,parent)
    % Write the <instancefields> section
    wrapper = document.createElement('InstanceFields');
    for i=1:size(self.xTrain,2)
        name = sprintf('x%i',i);
        instancefield = document.createElement('InstanceField');
        instancefield.setAttribute('field',name);
        instancefield.setAttribute('column',name);
        wrapper.appendChild(instancefield);
    end
    for i=1:size(self.yTrain,2)
        name = sprintf('y%i',i);
        instancefield = document.createElement('InstanceField');
        instancefield.setAttribute('field',name);
        instancefield.setAttribute('column',name);
        wrapper.appendChild(instancefield);
    end
    parent.appendChild(wrapper);
end

function writeInlineTable(self,document,parent)
    % Write <inlinetable> section that contains the training data
    nRows = size(self.xTrain,1);
    wrapper = document.createElement('InlineTable');
    for row=1:nRows
        rowElement = document.createElement('row');
        for i=1:size(self.xTrain,2)
            % Write the x columns
            name = sprintf('x%i',i);
            value = self.xTrain(row,i);
            node = document.createElement(name);
            node.setTextContent(num2str(value));
            rowElement.appendChild(node);
        end
        for i=1:size(self.yTrain,2)
            % Write the y columns
            name = sprintf('y%i',i);
            value = self.yTrain(row,i);
            node = document.createElement(name);
            node.setTextContent(num2str(value));
            rowElement.appendChild(node);
        end
        wrapper.appendChild(rowElement);
    end
    parent.appendChild(wrapper);
end





