<?xml version="1.0" ?>
<sld:StyledLayerDescriptor version="1.0.0" xmlns="http://www.opengis.net/sld" xmlns:gml="http://www.opengis.net/gml" xmlns:ogc="http://www.opengis.net/ogc" xmlns:sld="http://www.opengis.net/sld">
    <sld:UserLayer>
        <sld:LayerFeatureConstraints>
            <sld:FeatureTypeConstraint/>
        </sld:LayerFeatureConstraints>
        <sld:UserStyle>
            <sld:Name>ACDWRB</sld:Name>
            <sld:Title/>
            <sld:FeatureTypeStyle>
                <sld:Name/>
                <sld:Rule>
                    <sld:RasterSymbolizer>
                        <sld:Geometry>
                            <ogc:PropertyName>Grade of a sub-soil being acid e.g. having a low pH and low BS (SoilGrids250m)</ogc:PropertyName>
                        </sld:Geometry>
                        <sld:Opacity>1</sld:Opacity>
                        <sld:ColorMap>
                            <sld:ColorMapEntry color="#fff5eb" label="Low or none acidity problems" opacity="1.0" quantity="0"/>
                            <sld:ColorMapEntry color="#fdd1a5" label="Slightly acid (pH&lt;6.6) and/or dystric soil properties" opacity="1.0" quantity="1"/>
                            <sld:ColorMapEntry color="#fd9243" label="Moderately acid (pH&lt;5.5) and/or dystric Podzols and/or Gleysols" opacity="1.0" quantity="2"/>
                            <sld:ColorMapEntry color="#de4f05" label="Strongly acid (pH&lt;5) and/or dystric Acrisols and/or Alisols" opacity="1.0" quantity="3"/>
                            <sld:ColorMapEntry color="#7f2704" label="Extremely acid (pH&lt;4.5) and/or Acrisols, Alisols and/or alumic properties" opacity="1.0" quantity="4"/>
                        </sld:ColorMap>
                    </sld:RasterSymbolizer>
                </sld:Rule>
            </sld:FeatureTypeStyle>
        </sld:UserStyle>
    </sld:UserLayer>
</sld:StyledLayerDescriptor>
