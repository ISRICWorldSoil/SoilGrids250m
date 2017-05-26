<StyledLayerDescriptor xmlns="http://www.opengis.net/sld" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xmlns:ogc="http://www.opengis.net/ogc" xmlns:gml="http://www.opengis.net/gml" xsi:schemaLocation="http://www.opengis.net/sld StyledLayerDescriptor.xsd" version="1.0.0">
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
              <ogc:PropertyName>Predicted probability of occurrence (0-100%) of R horizon</ogc:PropertyName>
            </Geometry>
                        <ColorMap>
                            <ColorMapEntry color="#ffffd4" label="0.000000" opacity="1.0" quantity="0"/>
                            <ColorMapEntry color="#fee9ac" label="12.142857" opacity="1.0" quantity="12.1429"/>
                            <ColorMapEntry color="#fecf7f" label="24.285714" opacity="1.0" quantity="24.2857"/>
                            <ColorMapEntry color="#feab45" label="36.428571" opacity="1.0" quantity="36.4286"/>
                            <ColorMapEntry color="#f38821" label="48.571429" opacity="1.0" quantity="48.5714"/>
                            <ColorMapEntry color="#de6711" label="60.714286" opacity="1.0" quantity="60.7143"/>
                            <ColorMapEntry color="#bd4c09" label="72.857143" opacity="1.0" quantity="72.8571"/>
                            <ColorMapEntry color="#993404" label="85.000000" opacity="1.0" quantity="85"/>
                            <ColorMapEntry color="#FFFFFF" quantity="255" label="NODATA" opacity="0.0"/>
                        </ColorMap>
                    </RasterSymbolizer>
                </Rule>
            </FeatureTypeStyle>
        </UserStyle>
    </UserLayer>
</StyledLayerDescriptor>
