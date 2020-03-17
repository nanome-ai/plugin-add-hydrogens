if [ "$(docker ps -aq -f name=nanome-add-hydrogens)" != "" ]; then
    # cleanup
    echo "removing exited container"
    docker rm -f nanome-add-hydrogens
fi

docker run -d \
--name nanome-add-hydrogens \
--restart unless-stopped \
-e ARGS="$*" \
nanome-add-hydrogens
