<?xml version="1.0" ?>
<StyledLayerDescriptor version="1.0.0" xmlns="http://www.opengis.net/sld" xmlns:gml="http://www.opengis.net/gml" xmlns:ogc="http://www.opengis.net/ogc" xmlns:sld="http://www.opengis.net/sld">
    <UserLayer>
        <LayerFeatureConstraints>
            <FeatureTypeConstraint/>
        </LayerFeatureConstraints>
        <UserStyle>
            <Name>BDRLOG</Name>
            <Title/>
            <FeatureTypeStyle>
                <Name/>
                <Rule>
                    <RasterSymbolizer>
                        <Geometry>
              <PropertyName>Predicted probability of occurrence (0-100%) of R horizon</PropertyName>
              <Opacity>1</Opacity>
            </Geometry>
                        <ColorMap>
                        
                            <ColorMapEntry color="#ffffd4" label="0.000000" opacity="0.7" quantity="0"/>
                            <ColorMapEntry color="#fee9ac" label="12.142857" opacity="0.7" quantity="12.1429"/>
                            <ColorMapEntry color="#fecf7f" label="24.285714" opacity="0.7" quantity="24.2857"/>
                            <ColorMapEntry color="#feab45" label="36.428571" opacity="0.7" quantity="36.4286"/>
                            <ColorMapEntry color="#f38821" label="48.571429" opacity="0.7" quantity="48.5714"/>
                            <ColorMapEntry color="#de6711" label="60.714286" opacity="0.7" quantity="60.7143"/>
                            <ColorMapEntry color="#bd4c09" label="72.857143" opacity="0.7" quantity="72.8571"/>
                            <ColorMapEntry color="#993404" label="85.000000" opacity="0.7" quantity="85"/>
                            <ColorMapEntry color="#FFFFFF" quantity="255" label="NODATA" opacity="0.0"/>
                        </ColorMap>
                    </RasterSymbolizer>
                </Rule>
            </FeatureTypeStyle>
        </UserStyle>
    </UserLayer>
</StyledLayerDescriptor>
