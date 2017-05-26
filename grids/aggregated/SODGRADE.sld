<sld:StyledLayerDescriptor xmlns="http://www.opengis.net/sld" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xmlns:ogc="http://www.opengis.net/ogc" xmlns:gml="http://www.opengis.net/gml" xsi:schemaLocation="http://www.opengis.net/sld StyledLayerDescriptor.xsd" version="1.0.0">
    <sld:UserLayer>
        <sld:LayerFeatureConstraints>
            <sld:FeatureTypeConstraint/>
        </sld:LayerFeatureConstraints>
        <sld:UserStyle>
            <sld:Name>SLGWRB</sld:Name>
            <sld:Title/>
            <sld:FeatureTypeStyle>
                <sld:Name/>
                <sld:Rule>
                    <sld:RasterSymbolizer>
                        <sld:Geometry>
                            <ogc:PropertyName>Sodic soil grade based on WRB soil types and soil pH (SoilGrids250m)</ogc:PropertyName>
                        </sld:Geometry>
                        <sld:Opacity>1</sld:Opacity>
                        <ColorMap>
                            <ColorMapEntry color="#f7f4f9" label="No sodic properties" opacity="1.0" quantity="0"/>
                            <ColorMapEntry color="#d4b9da" label="Some sodic properties" opacity="1.0" quantity="1"/>
                            <ColorMapEntry color="#df65b0" label="Moderate sodic properties" opacity="1.0" quantity="2"/>
                            <ColorMapEntry color="#ce1256" label="Solonetz or Calcisols, sodic properties and/or ph &gt;8.1" opacity="1.0" quantity="3"/>
                            <ColorMapEntry color="#67001f" label="Solonetz, Solonchak, gypsic properties and/or ph &gt;8.5" opacity="1.0" quantity="4"/>
                            <ColorMapEntry color="#FFFFFF" quantity="255" label="NODATA" opacity="0.0"/>
                        </ColorMap>
                    </sld:RasterSymbolizer>
                </sld:Rule>
            </sld:FeatureTypeStyle>
        </sld:UserStyle>
    </sld:UserLayer>
</sld:StyledLayerDescriptor>
