'''
Created on Dec 22, 2011

@author: aaron
'''
import subprocess
import urllib2
import re
import os

class FirehoseConnection:
    firehose_base = "http://firehose:8080/cga/ws/"
    list_base = "list/"
    annot_base = "entity/getAnnotations/"
    sample_entity = "Sample?"
    individual_entity = "Individual?"
    filter_set_type = "filterSetType={filterSet}"
    filter_set_name = "filterSetName={filterName}"
    workspace_name = "workspaceName={workspace}"
    get_container = "getContainers="
    sep = "&"
    container_type = "containerType={container}"
    entity_name = "entityNames={entry}"

    def __init__(self,username,password):
        # save off the username and password
        self.username = username
        self.password = password

        # setup a password manager
        password_mgr = urllib2.HTTPPasswordMgrWithDefaultRealm()

        # Add the username and password.
        # If we knew the realm, we could use it instead of ``None``.
        password_mgr.add_password(None, FirehoseConnection.firehose_base, username, password)

        handler = urllib2.HTTPBasicAuthHandler(password_mgr)

        # create "opener" (OpenerDirector instance)
        opener = urllib2.build_opener(handler)

        # Install the opener.
        urllib2.install_opener(opener)

    def createCommand(self,workspace,individual_set,sample=True,get_container=True,sample_name=""):
        command = FirehoseConnection.firehose_base + FirehoseConnection.list_base
        if sample:
            command += FirehoseConnection.sample_entity
        else:
            command += FirehoseConnection.individual_entity

        if get_container:
            command += FirehoseConnection.get_container + FirehoseConnection.sep + FirehoseConnection.container_type.format(container="Individual") + FirehoseConnection.sep
            command += FirehoseConnection.workspace_name.format(workspace=workspace) + FirehoseConnection.sep + FirehoseConnection.entity_name.format(entry=sample_name)
        else:
            command += FirehoseConnection.filter_set_type.format(filterSet="Individual_Set") + FirehoseConnection.sep + FirehoseConnection.filter_set_name.format(filterName=individual_set)
            command += FirehoseConnection.sep + FirehoseConnection.workspace_name.format(workspace=workspace)
        return command

    def getSampleList(self,workspace,individual_set):
        command = self.createCommand(workspace,individual_set,get_container=False)
        resp = urllib2.urlopen(command)
        html = resp.read().split("\n")
        samples = [ln.strip() for ln in html]
        return(samples)

    def getSampleBams(self,workspace,individual_set,sample):
        command = FirehoseConnection.firehose_base + FirehoseConnection.annot_base + FirehoseConnection.sample_entity
        command += FirehoseConnection.entity_name.format(entry=sample) # + FirehoseConnection.sep # + FirehoseConnection.filter_set_type.format(filterSet="Individual_Set")
        command += FirehoseConnection.sep + FirehoseConnection.workspace_name.format(workspace=workspace) + FirehoseConnection.sep
        command += "annotationTypes=clean_bam_file_capture" + FirehoseConnection.sep + "annotationTypes=sample_type"
        resp = urllib2.urlopen(command)
        html = resp.read().split("\n")
        return(html[1].split("\t")[1:len(html[1].split("\t"))])

    def getIndividualList(self,workspace,individual_set):
        command = self.createCommand(workspace,individual_set,get_container=False,sample=False)
        resp = urllib2.urlopen(command)
        html = resp.read().split("\n")
        indv = [ln.strip() for ln in html]
        return(indv)

    def getMapping(self,workspace,individual_set,sample,individuals):
        command = self.createCommand(workspace,individual_set,get_container=True,sample_name=sample)
        resp = urllib2.urlopen(command)
        html = resp.read().split("\n")
        if not html[0].startswith("#Sample: "):
            print "zero = " + html[0]
            print html
            raise NameError("Unable to parse command " + command)

        ind = []
        for i in range(1,len(html)):
            if html[i] != "":
                if html[i] in individuals:
                    return (html[i],sample)
        return (None,None)

    def getIndividualSampleMapping(self,workspace,individual_set):
        indvs = self.getIndividualList(workspace,individual_set)
        samples = self.getSampleList(workspace,individual_set)
        mapping = []
        for smp in filter(lambda x: x != "",samples):
            mp = self.getMapping(workspace,individual_set,smp,indvs)
            bam_type = self.getSampleBams(workspace,individual_set,smp)
            if (mp[0] != None):
                mapping.append((mp[0],mp[1],bam_type[0],bam_type[1]))
        return(mapping)

    def resolve_picard_path(self,picard_path):
        if not picard_path.startswith("/seq/picard_aggregation/"):
            if not os.path.exists(picard_path):
                return "NA"
            else:
                return picard_path
        else:
            path_split = [re.sub('v\d+', 'current',token) for token in picard_path.split("/")]
            pth = "/".join(path_split)
            if not os.path.exists(pth):
                return "NA"
            else:
                return pth

    def mappingToOutputFile(self,workspace,individual_sets,output_file):
        of = open(output_file,"w")
        mapping = []
        for ind_set in individual_sets:
            mapping.extend(self.getIndividualSampleMapping(workspace,ind_set))
        individuals = set()
        for it in mapping:
            individuals.add(it[0])
        nbams = set()
        tbams = set()
        of.write("sample\ttumor_bam\tnormal_bam\n")
        for ind in individuals:
            tumor_bam = "NA"
            normal_bam = "NA"
            for ln in mapping:
                if ln[0] == ind and ln[3] == "Tumor":
                    tumor_bam = ln[2]
            for ln in mapping:
                if ln[0] == ind and ln[3] == "Normal":
                    normal_bam = ln[2]

            if tumor_bam != "NA" and not tumor_bam in tbams:
                tbams.add(tumor_bam)
            else:
                tumor_bam = "NA"
            if normal_bam != "NA" and not normal_bam in nbams:
                nbams.add(normal_bam)
            else:
                normal_bam = "NA"
            of.write(ind + "\t" + self.resolve_picard_path(tumor_bam) + "\t" + self.resolve_picard_path(normal_bam) + "\n")
