default:
	@echo pull chr19

pull:
	docker pull sagebionetworks/synapsepythonclient

bash:
	docker run -it \
           --entrypoint /bin/bash \
           -v /tmp/docker-paul:/tmp \
           sagebionetworks/synapsepythonclient

khaleesi-bash:
	sudo docker run -it \
           --entrypoint /bin/bash \
           -v /tmp/docker-paul:/tmp \
           sagebionetworks/synapsepythonclient


khaleesi-bash-new:
	sudo docker run -it \
           --entrypoint /bin/bash \
           -v /tmp/docker-paul:/tmp \
	    docker.synapse.org/syn25326461/synpy-1128


#------------- notes
# export PS1='DOCKER.\W> '
# python3
# import synapseclient
# syn = synapseclient.login("paul-shannon", password=oldAnimal)
# x = syn.get('syn1899498', downloadLocation="/tmp")  # quick and easy
# x = syn.get('syn11714133', downloadLocation="/tmp")


