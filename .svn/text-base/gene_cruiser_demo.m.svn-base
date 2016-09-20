cd ~/matlab

% Create class using the wsdl
wsdlURL = 'http://www.broad.mit.edu/webservices/genecruiser/services/Annotation?wsdl';
className = createClassFromWsdl(wsdlURL);
% This creates a @AnnotationServiceService directory under the current directory
% with functions for the specific class

% Check the available methods
methods(AnnotationServiceService);

% create an instance of the class
ass=AnnotationServiceService;
% check the class is indeed correct
class(ass)

dbmap=getDatabaseToFieldsMap(ass);
display_soap_result(dbmap);

help annotateProbes

% Example from Bioinformatics Application Note
query_fields={'LocusLink Locus ID','Unigene_Human UniGene Cluster'};
probes={'urn:lsid:affymetrix.com:probeset.hu6800:AF000430_at','AB002409_at'};
annot_res=annotateProbes(ass,probes,query_fields);
